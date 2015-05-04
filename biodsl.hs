import System.Environment
import Control.Monad
import Data.List
import Data.Algorithms.KMP
import Alignment
import Text.PrettyPrint
import Text.PrettyPrint.HughesPJClass

-- Types of data to operate on
data Nucleotide = A | C | G | T | N deriving (Eq,Show)
type DNA = [Nucleotide]

data Feature = Nil | Single (Int,Int) | Multipart (Int,Int) [Feature]  deriving (Eq,Show)

instance Pretty Feature where
	pPrint Nil = text "Nil"
	pPrint (Single (a,b)) = text "Single (" <> pPrint a <> text "," <> pPrint b <> text ")"
	pPrint (Multipart (a,b) parts) = vcat $ [ text "Multipart (" <> pPrint a <> text "," <> pPrint b <> text ")" ] ++ map (nest 3 . pPrint) parts

data Signal
	= Signal (DNA -> Feature)
	| Signal :> Signal
	| (Signal,Int) :#: Signal

motif :: String -> Signal
motif m =
  let kmpTable = build $ dna m
  in Signal (\d -> 
				let matches = match kmpTable d
				in if (length matches) > 0 then Single ((head matches)+1,(head matches)+(length m)) else Nil)

invertedRepeats :: Signal
invertedRepeats = Signal (\d ->
  	let
      aligns = (align cost d (reverseComplement d))
      l = length d
		in if (length aligns) == 0 
				then Nil
				else let (p,q,r,s) = (head aligns)
							in Multipart (p,l-r+1) [Single (1,q-p+1), Single (q-p+1,l-s-p+2), Single (l-s-p+2,l-r-p+2)])

gcRich :: Signal
gcRich = Signal (\d -> if (toRational (length [x | x <- d, x == G || x == C])) > (0.75 * (toRational (length d))) then Single (1, (length d)) else Nil)

part :: Int -> Feature -> Feature
part _ Nil = Nil
part _ (Single _) = Nil
part i (Multipart (a,b) parts) = if (i >= 0 && i < (length parts)) then (scale a (parts!!i)) else Nil

slice from to xs = take (to - from + 1) (drop (from-1) xs)

replaceAtIndex n item ls = a ++ (item:b) where (a, (_:b)) = splitAt n ls

scale :: Int -> Feature -> Feature
scale i Nil = Nil
scale i (Single (a,b)) = Single (a+i-1,b+i-1)
scale i (Multipart (a,b) parts) = Multipart (a+i-1,b+i-1) parts

interval :: Feature -> Maybe (Int,Int)
interval Nil = Nothing
interval (Single inter) = Just inter
interval (Multipart inter parts) = Just inter

findFeature :: Signal -> DNA -> Feature
findFeature (Signal s) d = (s d) 
findFeature (a :> b) d
	= case (findFeature a d) of 
			Nil -> Nil
			Single (start1,end1) -> case (findFeature b (drop end1 d)) of
															Nil -> Nil
															Single (start2,end2) -> Multipart (start1,end1+end2) 
																												[Single (1,end1-start1+1), 
																												Single (end1-start1+1,end1+start2-start1+1), 
																												Single (end1+start2-start1+1,end1+end2-start1+1)]
															Multipart (start2,end2) parts -> Multipart (start1,end1+end2) 
																																[Single (1,end1-start1+1), 
																																Single (end1-start1+1,end1+start2-start1+1), 
																																Multipart (end1+start2-start1+1,end1+end2-start1+1) parts]
			Multipart (start1,end1) parts1 -> case (findFeature b (drop end1 d)) of
                              Nil -> Nil
                              Single (start2,end2) -> Multipart (start1,end1+end2) 
																												[Multipart (1,end1-start1+1) parts1, 
																												Single (end1-start1+1,end1+start2-start1+1), 
																												Single (end1+start2-start1+1,end1+end2-start1+1)]
                              Multipart (start2,end2) parts2 -> Multipart (start1,end1+end2) 
																																	[Multipart (1,end1-start1+1) parts1, 
																																	Single (end1-start1+1,end1+start2-start1+1), 
																																	Multipart (end1+start2-start1+1,end1+end2-start1+1) parts2]

findFeature ((a,i) :#: b) d
	= case (findFeature a d) of
			Nil -> Nil
			Single (start1,end1) -> if i /= 0 
															then Nil 
															else case (findFeature b (slice start1 end1 d)) of
																Nil -> Nil
																Single (start2,end2) -> Multipart (start1,end1) 
																													[Single (1,start2), 
																													Single (start2,end2), 
																													Single (end2,end1-start1+1)]
																Multipart (start2,end2) parts -> Multipart (start1,end1)
																																	[Single (1,start2),
																																	Multipart (start2,end2) parts,
																																	Single (end2,end1-start1+1)]
			Multipart (start1,end1) parts1 -> 
				case interval (part i (Multipart (start1,end1) parts1)) of
					Nothing -> Nil
					Just (start2,end2) -> case (findFeature b (slice start2 end2 d)) of
															Nil -> Nil
															Single (start3,end3) -> Multipart (start1,end1) 
																												(replaceAtIndex i (Multipart (start2-start1+1,end2-start1+1)
																																						[Single (1,start3),
																																						Single (start3,end3),
																																						Single (end3,end2-start2+1)]) parts1)
															Multipart (start3,end3) parts -> Multipart (start1,end1)
																																(replaceAtIndex i (Multipart (start2-start1+1,end2-start1+1)
																																										[Single (1,start3),
																																										Multipart (start3,end3) parts,
																																										Single (end3,end2-start2+1)]) parts1)
																


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
	putStrLn $ prettyShow $ findFeature invertedRepeats (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ findFeature ((invertedRepeats,0) :#: gcRich) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ findFeature ((((invertedRepeats,0) :#: gcRich),1) :#: transposase) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	putStrLn $ prettyShow $ findFeature ((((invertedRepeats,0) :#: gcRich),1) :#: (transposase :> pif)) (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAGGACGGCGGCGGCGTTATCTAGGTATCGACTCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")

