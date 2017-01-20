module Octell (
    -- * Matrix type
      Matrix

    -- * Query
    , numRows
    , numCols
    , size
    , isScalar
    , isRowVector
    , isColVector
    , isSquare
    , get

    -- * Construction
    , identity
    , gen
    , zeros
    , allN
    , fromList
    , randn

    -- ** Deconstruction
    , toList

    -- ** Manipulation
    , row
    , col
    , dropRow
    , dropCol
    , unionRows
    , unionCols

    -- * Matrix operations
    , det
    , inv
    , transpose
    , trace
    , minor
    , cofactor
    ) where

import Control.Monad (replicateM)
import Data.List (groupBy)
import Data.Map (Map)
import qualified Data.Map as M
import System.Random (randomIO, Random)
import System.IO.Unsafe (unsafePerformIO)

data Matrix e = Matrix
    { elements :: Map (Int, Int) e
    } deriving (Eq)

instance Show e => Show (Matrix e) where
    show m
      | isScalar m = "scalar" ++ unlines show'
      | isRowVector m = "row vector (" ++ show (numCols m) ++ ")" ++ unlines show'
      | isColVector m = "column vector (" ++ show (numRows m) ++ ")" ++ unlines show'
      | otherwise = "matrix (" ++ show (numRows m) ++ " x " ++ show (numCols m) ++ ")" ++ unlines show'
      where mat = map (map show) (toList m)
            longest = maximum . map length $ concat mat
            pad x = x ++ replicate (longest - length x) ' '
            show' = "" : "" : map ((:) '\t' . unwords . map pad) mat

instance Functor Matrix where
    fmap f m = m { elements = M.map f (elements m) }

instance (Num e, Show e) => Num (Matrix e) where
    m + m' = byElem (+) m m'
    m * m' = prod m m'
    abs = fmap abs
    negate = fmap negate
    signum = fmap signum
    fromInteger i = fromList [[fromInteger i]]

instance (Eq e, Show e, Fractional e) => Fractional (Matrix e) where
    recip = inv
    fromRational r = fromList [[fromRational r]]

instance (Show e, Eq e, Floating e) => Floating (Matrix e) where
    pi = fromList [[pi]]
    exp = fmap exp
    log = fmap log
    sin = fmap sin
    cos = fmap cos
    asin = fmap asin
    acos = fmap acos
    atan = fmap atan
    sinh = fmap sinh
    cosh = fmap cosh
    asinh = fmap asinh
    acosh = fmap acosh
    atanh = fmap atanh

-- | Visualise a matrix in a very rudimentary way
plotMatrix :: (Fractional e, Num e, Ord e) => Matrix e -> IO ()
plotMatrix m = putStrLn . unlines . map concat . toList $ f <$> m
    where f x | x < minVal + range   = "  "
              | x < minVal + 2*range = "░░"
              | x < minVal + 3*range = "▒▒"
              | x < minVal + 4*range = "▓▓"
              | otherwise            = "██"
          maxVal = maximum . concat $ toList m
          minVal = minimum . concat $ toList m
          range = (maxVal - minVal)/5

-- | Get the number of rows of a given matrix
numRows :: Matrix e -> Int
numRows (Matrix m) = maximum . map fst $ M.keys m

-- | Get the number of columns of a given matrix
numCols :: Matrix e -> Int
numCols (Matrix m) = maximum . map snd $ M.keys m

-- | Get the (i,j) size of a given matrix
size :: Matrix e -> (Int, Int)
size m = (numRows m, numCols m)

-- | Determine if a matrix is just a scalar value
isScalar :: Matrix e -> Bool
isScalar (Matrix m) = M.size m == 1

-- | Determine if a matrix is just a row vector
isRowVector :: Matrix e -> Bool
isRowVector m = numRows m == 1

-- | Determine if a matrix is just a column vector
isColVector :: Matrix e -> Bool
isColVector m = numCols m == 1

-- | Determine if a matrix is square
isSquare :: Matrix e -> Bool
isSquare m = numRows m == numCols m

-- | Make a matrix from a list of lists
fromList :: [[e]] -> Matrix e
fromList es = Matrix . M.fromList $ zip indices trimmed
    where minCols = minimum $ map length es
          trimmed = concatMap (take minCols) es
          indices = [ (x,y) | x <- [1..], y <- [1..minCols] ]

-- | Convert matrix to list of lists
toList :: Matrix e -> [[e]]
toList = map (map snd) . groupBy f . M.toAscList . elements
    where f x y = fst (fst x) == fst (fst y)

-- | Make an i x j matrix of given value
allN :: (Int, Int) -> e -> Matrix e
allN (i,j) = fromList . replicate i . replicate j

-- | Make an i x j matrix of zeros
zeros :: Num e => (Int, Int) -> Matrix e
zeros (i,j) = allN (i,j) 0

-- | Make an i x j matrix using the given generator function
gen :: (Int, Int) -> ((Int, Int) -> e) -> Matrix e
gen (i,j) f = Matrix . M.mapWithKey g . elements $ zeros (i,i)
    where g a _ = f a

-- | Make an i x i identity matrix
identity :: Num e => Int -> Matrix e
identity i = gen (i,i) f
    where f (x,y) = if x == y then 1 else 0

-- | Naughty random matrix
randn :: (Int, Int) -> Matrix Double
randn (i,j) = fromList . unsafePerformIO . replicateM i $ replicateM j randomIO

-- | Perform an element wise binary operation on two matrices
byElem :: Show e => (e -> e -> e) -> Matrix e -> Matrix e -> Matrix e
byElem f m m' = if numCols m /= numCols m' || numRows m /= numRows m'
                   then error $ "Matrices have different dimensions when "
                                ++ "attempting binary elementwise operation\n"
                                ++ show m ++ show m'
                   else Matrix $ M.unionWith f (elements m) (elements m')

-- | Transpose a matrix
transpose :: Matrix e -> Matrix e
transpose m = Matrix $ M.mapKeys (\(x,y) -> (y,x)) (elements m)

-- | Calculate the trace of a matrix
trace :: (Show e, Num e) => Matrix e -> e
trace m
  | not (isSquare m) = error $ "Can't take the trace of non square"
                             ++ " matrix:\n" ++ show m
  | otherwise = sum . M.filterWithKey (\(x,y) _ -> x == y) $ elements m

-- | Calculate the determinant of a matrix
det :: (Show e, Num e) => Matrix e -> e
det m
  | not (isSquare m) = error $ "Can't take the determinant of non square"
                             ++ " matrix:\n" ++ show m
  | numRows m < 2 = error $ "Determinant requires at least a 2x2 matrix:\n"
                          ++ show m
  | numRows m == 2 = (get (1,1) m) * (get (2,2) m)
                   - (get (1,2) m) * (get (2,1) m)
  | otherwise = sum $ zipWith (*) firstRow subDets
  where firstRow = zipWith (*) (head (toList m)) (cycle [1, -1])
        indices = [ (1,x) | x <- [1..] ] -- list of first row indices
        subDets = map (flip minor m) indices -- list of subdeterminants

-- | Naive matrix product
prod :: (Show e, Num e) => Matrix e -> Matrix e -> Matrix e
prod m m'
  | isScalar m = (* get (1,1) m) <$> m'
  | isScalar m' = (* get (1,1) m') <$> m
  | numCols m /= numRows m' = error $ "Can't take matrix product of:\n"
                                    ++ unlines [show m, show m']
  | otherwise = fromList $ splitEvery (numCols m') apList
  where apList = map sum $ (zipWith (*))
              <$> toList m
              <*> toList (transpose m')

-- | Helper function, split list every i items
splitEvery :: Int -> [a] -> [[a]]
splitEvery i [] = [[]]
splitEvery i xs = if length xs <= i
                    then [take i xs]
                    else take i xs : splitEvery i (drop i xs)

-- | Matrix inverse by Cramer's rule
inv :: (Eq e, Num e, Show e, Fractional e) => Matrix e -> Matrix e
inv m
  | isScalar m = recip <$> m
  | not (isSquare m) = error $ "Can't take the inverse of non square"
                             ++ " matrix:\n" ++ show m
  | det m == 0 = error $ "Determinant is 0, inverse does not exist:\n"
                       ++ show m
  | otherwise = (/ (det m)) <$> transpose (fromList cofactors)
  where cofactors = map (map (flip cofactor m)) indices
        indices = [ [ (r,c) | c <- [1..numCols m] ] | r <- [1..numRows m] ]

-- | Get the (i,j)th element of a matrix
get :: Show e => (Int, Int) -> Matrix e -> e
get (i,j) m
  | i <= 0 || j <= 0 = error $ "Row and column must be positive integers, "
                             ++ "given " ++ show (i,j)
  | otherwise = case M.lookup (i,j) $ elements m of
                  Just e -> e
                  Nothing -> error $ "Specified row and column " ++ show (i,j)
                                   ++ " are out of bounds:\n" ++ show m

-- | First minor of matrix (only deletes one row and column).
minor :: (Num e, Show e) => (Int, Int) -> Matrix e -> e
minor (i,j)
  | i <= 0 || j <= 0 = error $ "Row and column must be positive integers, "
                             ++ "given " ++ show (i,j)
  | otherwise = det . dropCol j . dropRow i

cofactor :: (Num e, Show e) => (Int, Int) -> Matrix e -> e
cofactor (i,j) m = (-1)^(i+j) * minor (i,j) m

-- | Drop a row from the matrix
dropRow :: Show e => Int -> Matrix e -> Matrix e
dropRow i m
  | i <= 0 = error $ "Row must be a positive integer, given " ++ show i
  | numRows m < i = m
  | otherwise = Matrix . M.mapKeys f . M.filterWithKey g $ elements m
  where f (r,c) = if r > i then (r-1,c) else (r,c)
        g (r,_) _ = r /= i

-- | Drop a column from the matrix
dropCol :: Show e => Int -> Matrix e -> Matrix e
dropCol j m
  | j <= 0 = error $ "Col must be a positive integer, given " ++ show j
  | numCols m < j = m
  | otherwise = Matrix . M.mapKeys f . M.filterWithKey g $ elements m
  where f (r,c) = if c > j then (r,c-1) else (r,c)
        g (_,c) _ = c /= j

-- | Make a row vector from a specified row of the matrix
row :: Show e => Int -> Matrix e -> Matrix e
row i m
  | i <= 0 = error $ "Row must be a positive integer, given " ++ show i
  | numRows m < i = error $ "Can't get row " ++ show i
                          ++ " of smaller matrix:\n" ++ show m
  | otherwise = Matrix . M.mapKeys f . M.filterWithKey g $ elements m
  where f (_,c) = (1,c)
        g (r,_) _ = r == i

-- | Make a column vector from a specified column of the matrix
col :: Show e => Int -> Matrix e -> Matrix e
col j m
  | j <= 0 = error $ "Col must be a positive integer, given " ++ show j
  | numCols m < j = error $ "Can't get col " ++ show j
                          ++ " of smaller matrix:\n" ++ show m
  | otherwise = Matrix . M.mapKeys f . M.filterWithKey g $ elements m
  where f (r,_) = (r,1)
        g (_,c) _ = c == j

-- | Add rows of second matrix onto rows of first matrix
unionRows :: Show e => Matrix e -> Matrix e -> Matrix e
unionRows m m'
  | numCols m /= numCols m' = error $ "Can't join matrices of different "
                                    ++ "column counts:\n"
                                    ++ unlines [show m, show m']
  | otherwise = Matrix . M.union (elements m) $ M.mapKeys f (elements m')
  where f (r,c) = (r + numRows m, c)

-- | Add columns of second matrix onto columns of first matrix
unionCols :: Show e => Matrix e -> Matrix e -> Matrix e
unionCols m m'
  | numRows m /= numRows m' = error $ "Can't join matrices of different "
                                    ++ "row counts:\n"
                                    ++ unlines [show m, show m']
  | otherwise = Matrix . M.union (elements m) $ M.mapKeys f (elements m')
  where f (r,c) = (r, c + numCols m)
