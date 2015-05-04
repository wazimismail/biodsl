module Alignment where
import Data.Array
import Data.List
import Data.Function (on)
import Data.Ord

data Action = Match | Insertion | Deletion | Mismatch | None deriving (Show,Eq)

type Distance = Int

tfst :: (a, b, c) -> a
tfst (a,_,_) = a

thd :: (a, b, c) -> c
thd (_,_,c) = c


align :: Eq a => (Action -> Distance) -> [a] -> [a] -> [(Int,Int,Int,Int)]
-- (Distance,[Action],[(Int,Int)]) 
-- ((Int,Int),(Int,Int))
align cost a b =
	let	
		(_,q,r) = (d m n)
		paths = ((filter ((== Match) . fst . head)) . (filter ((> 5) . length)) . (groupBy ((==) `on` fst))) (zip (reverse q) (reverse r))
	in 
		map (\path	-> ((fst $ snd $ head $ path)+1,(fst $ snd $ last $ path)+1,(snd $ snd $ head $ path)+1,(snd $ snd $ last $ path)+1))
				paths
  where (m, n) = (length a, length b)
        a'     = listArray (1, m) a
        b'     = listArray (1, n) b

        d 0 0 = (0, [], [])
        d i 0 = go (i - 1) 0 Insertion
        d 0 j = go 0 (j - 1) Deletion
        d i j
          | a' ! i ==  b' ! j = go (i - 1) (j - 1) Match
          | otherwise = maximum' [ go (i - 1) j       Deletion
                                 , go i (j - 1)       Insertion
                                 , go (i - 1) (j - 1) Mismatch
																 , go (i - 1) (j - 1) None
                                 ]

        maximum' = maximumBy (comparing tfst)
        go i j action = let (score, actions, cells) = ds ! (i, j) in
          (score + cost action, (action : actions), ((i,j) : cells))

        ds = listArray bounds [d i j | (i, j) <- range bounds]
        bounds = ((0, 0), (m, n))

cost :: Action -> Distance
cost Match 		= 10
cost Mismatch = -5
cost None			= 0
cost _    		= -4
