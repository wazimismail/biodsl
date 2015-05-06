import System.Environment
import Control.Monad.Identity
import Control.Monad.Error
import Control.Monad.State
import Data.Maybe
import Data.List
import Data.Algorithms.KMP
import Alignment
import Text.PrettyPrint
import Text.PrettyPrint.HughesPJClass

-- Types of data to operate on
data Nucleotide = A | C | G | T | N deriving (Eq,Show)
type DNA = [Nucleotide]

data Feature = Node (Int,Int) [Feature]  deriving (Eq,Show)

instance Pretty Feature where
	pPrint (Node (a,b) parts) = vcat $ [ text "Node (" <> pPrint a <> text "," <> pPrint b <> text ")" ] ++ map (nest 3 . pPrint) parts

type FindFeature = (ErrorT String Identity) Feature

data Signal
	= Signal (DNA -> Maybe Feature)
	| Signal :> Signal
	| Part Int Signal
	| (Signal,Int) :#: Signal

motif :: String -> Signal
motif m =
  let kmpTable = build $ dna m
  in Signal (\d -> 
				let matches = match kmpTable d
				in if (length matches) > 0 then Just (Node ((head matches)+1,(head matches)+(length m)) []) else Nothing)

invertedRepeats :: Signal
invertedRepeats = Signal (\d ->
  	let
      aligns = (align cost d (reverseComplement d))
      l = length d
		in if (length aligns) == 0 
				then Nothing
				else let (p,q,r,s) = (head aligns)
							in Just $ Node (p,l-r+1) [Node (1,q-p+1) [], Node (q-p+1,l-s-p+2) [], Node (l-s-p+2,l-r-p+2) []])

gcRich :: Signal
gcRich = Signal (\d -> if (toRational (length [x | x <- d, x == G || x == C])) > (0.75 * (toRational (length d))) then Just $ Node (1, (length d)) [] else Nothing)

slice from to xs = take (to - from + 1) (drop (from-1) xs)

replaceAtIndex n item ls = a ++ (item:b) where (a, (_:b)) = splitAt n ls

runFindFeature :: FindFeature -> (Either String Feature)
runFindFeature ff = runIdentity (runErrorT ff)

findFeature :: Signal -> DNA -> FindFeature
findFeature (Signal s) d = case s d of
														Nothing -> throwError ("No feature found!")
														Just val -> return val

findFeature (a :> b) d = do
	(Node (start1,end1) parts1) <- (findFeature a d)
	(Node (start2,end2) parts2) <- (findFeature b (drop end1 d))
	return $ Node (start1,end1+end2) [Node (1,end1-start1+1) parts1, 
																		Node (end1-start1+1,end1+start2-start1+1) [], 
																		Node (end1+start2-start1+1,end1+end2-start1+1) parts2]

findFeature (Part i s) d = do
	(Node (start,end) parts) <- (findFeature s d)
	if (i < 0 || i >= (length parts)) 
	then throwError ("No feature found!")
	else 
		let (Node (starti,endi) partsi) = parts!!i in
		return $ (Node (starti+start-1,endi+start-1) partsi)

findFeature ((a,i) :#: b) d = do
	(Node (start1,end1) parts1) <- (findFeature a d)
	if (parts1 == []) 
	then do
		(Node (start2,end2) parts2) <- (findFeature b (slice start1 end1 d))
		return $ Node (start1,end1) [Node (1,start2) [], Node (start2,end2) parts2, Node (end2,end1-start1+1) []]
	else if (i < 0 || i >= (length parts1))
	then do
		throwError ("No feature found!")
	else
		let 
			(Node (starti,endi) partsi) = parts1!!i
			start2 = starti+start1-1
			end2 = endi+start1-1
		in do
			(Node (start3,end3) parts3) <- (findFeature b (slice start2 end2 d))
			return $ Node (start1,end1) (replaceAtIndex i (Node (start2-start1+1,end2-start1+1) 
																											[Node (1,start3) [], 
																											Node (start3,end3) parts3, 
																											Node (end3,end2-start2+1) []]) parts1)														


-- Helper functions to convert Char/String to Nucleotide/DNA
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


-- Test code
main :: IO ()
main = do 
	args <- getArgs
	-- let x = dna (args !! 0)
	let transposase = (motif "AGGTATCGACTCGTTTACGG")
	let pif = (motif "GAGCTGGGCTTT")
	let fwd = (dna "GGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCG")
	let rev = (dna "CGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCC")
	putStrLn $ show $ (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
--	putStrLn $ prettyShow $ runFindFeature $ findFeature (motif "ATG") (dna "TTCTTAGATGCGGC")
--	putStrLn $ prettyShow $ runFindFeature $ findFeature ((motif "ATG") :> (motif "GGC")) (dna "TTCTTAGATGCGGC")
	putStrLn $ prettyShow $ runFindFeature $ findFeature invertedRepeats (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
--	putStrLn $ prettyShow $ runFindFeature $ findFeature (Part 0 invertedRepeats) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ runFindFeature $ findFeature ((invertedRepeats,0) :#: gcRich) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ runFindFeature $ findFeature ((((invertedRepeats,0) :#: gcRich),1) :#: transposase) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ runFindFeature $ findFeature ((((invertedRepeats,0) :#: gcRich),1) :#: (transposase :> pif)) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")

