module Biodsl where
import Control.Monad
import Text.PrettyPrint
import Text.PrettyPrint.HughesPJClass

-- DSL
data Nucleotide = A | C | G | T | N deriving (Eq,Show)
type DNA = [Nucleotide]
data Feature = Node String (Int,Int) [Feature]  deriving (Eq,Show)
type Signal = DNA -> [Feature]
type State = (Int,Int)

instance Pretty Feature where
	pPrint (Node s (a,b) parts) = vcat $ [ pPrint s <> text " (" <> pPrint a <> text "," <> pPrint b <> text ")" ] ++ map (nest 3 . pPrint) parts

-- Monad definition
newtype SignalM a = SignalM (DNA -> State -> [(a,State)])
search (SignalM s) = s

liftF :: Signal -> SignalM Feature
liftF s = SignalM (\d (from,to) -> 
					(map 
						(\(Node s (a,b) c) -> ((Node s (a+from-1,b+from-1) c), (b+from-1,to))) 
						(s (slice d from to))))
				
instance Monad SignalM where
	return a	= SignalM (\d t -> [(a,t)])
	s >>= f		= SignalM (\d t1 -> concat [search (f a) d t2 | (a,t2) <- search s d t1])
	
put :: State -> SignalM ()
put newState = SignalM (\d t -> [((),newState)])

getValue :: (Feature, State) -> Feature
getValue (val,s) = val

getValues :: [(Feature, State)] -> [Feature]
getValues = map getValue
									
instance MonadPlus SignalM where
	mzero = SignalM (\d t -> [])
	p `mplus` q = SignalM (\d t -> search p d t ++ search q d t)

-- Combinators
dplus :: (Eq a) => SignalM a -> SignalM a -> SignalM a
p `dplus` q = SignalM (\d t -> case (search (p `mplus` q) d t) of
								[] -> []
								(c:cs) -> [c])

bplus :: (Eq a) => SignalM a -> SignalM a -> SignalM a
p `bplus` q = SignalM (\d t -> if ((search p d t) /= []) then (search p d t) else (search q d t))

sat :: SignalM a -> (a -> Bool) -> SignalM a
sat s p = do {c <- s; if (p c) then (return c) else mzero}

followedBy :: SignalM Feature -> SignalM Feature -> SignalM Feature
p `followedBy` q = do {c <- p; d <- q; return (merge (c:(d:[])))}

many :: SignalM a -> SignalM [a]
many s = (many1 s) `mplus` (return [])

many1 :: SignalM a -> SignalM [a]
many1 s = (do {c <- s; d <- (many s); return (c:d)})

tandem :: SignalM Feature -> SignalM [Feature]
tandem s = (tandem1 s) `dplus` (return [])

tandem1 :: SignalM Feature -> SignalM [Feature]
tandem1 s = (do {
				(Node s1 (start1,end1) parts1) <- s; 
				d <- (tandem s); 
				case d of 
					[] -> return [(Node s1 (start1,end1) parts1)]
					((Node s2 (start2,end2) parts2):cs) -> 
						if start2 == end1+1 then 
							return ((Node s1 (start1,end1) parts1):((Node s2 (start2,end2) parts2):cs)) 
						else return []})

mergeM :: SignalM [Feature] -> SignalM Feature
mergeM s = (do {c <- s; return (merge c)})

mergeMany :: SignalM Feature -> SignalM Feature
mergeMany s = mergeM $ many $ s

mergeMany1 :: SignalM Feature -> SignalM Feature
mergeMany1 s = mergeM $ many1 $ s

part :: Int -> SignalM Feature -> SignalM Feature
part j s = do {
			(Node name (start,end) parts) <- s;
			let 
				(Node sj (startj,endj) partsj) = parts!!(j-1)
			in do return $ (Node (name++"_"++sj) (startj+start-1,endj+start-1) partsj)}

havingA :: SignalM Feature -> SignalM Feature -> SignalM Feature
a `havingA` b = do {
					(Node s1 (start1,end1) parts1) <- a;
					put $ (start1,end1);
					(Node s2 (start2,end2) parts2) <- b;
					return $ (Node s1 (start1,end1) [Node "leftFlank" (1,start2-start1+1) [], 
														Node s2 (start2-start1+1,end2-start1+1) parts2, 
														Node "rightFlank" (end2-start1+1,end1-start1+1) []])}
			
havingB :: (SignalM Feature,Int) -> SignalM Feature -> SignalM Feature
(a,j) `havingB` b = do {
					(Node s1 (start1,end1) parts1) <- a;
					if j > (length parts1) || j <= 0 then mzero else
					let 
						(Node sj (startj,endj) partsj) = parts1!!(j-1)
						start2 = startj+start1-1
						end2 = endj+start1-1
					in do
						put $ (start2,end2);
						(Node s3 (start3,end3) parts3) <- b;
						put $ (start1,end1);
						return $ (Node s1 (start1,end1) (replaceAtIndex (j-1) (Node sj (startj,endj) 
									[Node "leftFlank" (1,start3-start2+1) [], 
									Node s3 (start3-start2+1,end3-start2+1) parts3, 
									Node "rightFlank" (end3-start2+1,end2-start2+1) []]) parts1))}
							
-- Helper functions

nucleotide :: Char -> Nucleotide
nucleotide 'A' = A
nucleotide 'C' = C
nucleotide 'G' = G
nucleotide 'T' = T
nucleotide _   = N

dna :: String -> DNA
dna = map nucleotide

complement :: Nucleotide -> Nucleotide
complement A = T
complement T = A
complement C = G
complement G = C
complement N = N

reverseComplement :: DNA -> DNA
reverseComplement d = reverse (map complement d)

slice :: [a] -> Int -> Int -> [a]
slice l from (-1) = (drop (from-1) l)
slice l from to = take (to - from + 1) (drop (from-1) l)

replaceAtIndex :: Int -> a -> [a] -> [a]
replaceAtIndex n item ls = c ++ (item:d) where (c, (_:d)) = splitAt n ls

merge :: [Feature] -> Feature
merge (c:[]) = c
merge ((Node s1 (start1,_) _):(Node s2 (_,end2) _):cs) = merge ((Node (s1++" "++s2) (start1,end2) []):cs)
