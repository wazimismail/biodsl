import System.Environment
import Control.Monad
import Biodsl
import Data.Algorithms.KMP
import Alignment
import Text.PrettyPrint
import Text.PrettyPrint.HughesPJClass


motif :: String -> String -> Signal
motif m name =
  let kmpTable = build $ dna m
  in (\d -> 
		let matches = match kmpTable d
		in if (length matches) > 0
			then let
					a = (head matches)+1
					b = a+(length m)-1
				in [(Node name (a,b) [])] 
			else [])
			
motifM :: String -> String -> SignalM Feature
motifM m name = liftF $ motif m name
				
allMotifs :: String -> String -> Signal
allMotifs m name =
  let kmpTable = build $ dna m
  in (\d -> 
		let matches = match kmpTable d
		in if (length matches) > 0 
			then 
				(map (\i -> (Node name (i+1,i+(length m)) [])) matches)
			else [])
			
allMotifsM :: String -> String -> SignalM Feature
allMotifsM m name = liftF $ allMotifs m name

directRepeats :: Signal
directRepeats = (\d -> [Node "dRep" (10,80) [Node "repL" (1,5) [], Node "between" (5,67) [], Node "repR" (67,71) []]])

directRepeatsM :: SignalM Feature
directRepeatsM = liftF $ directRepeats

invertedRepeats :: Signal
invertedRepeats = (\d -> ((Node "iRep" (2,62) [Node "repL" (1,10) [], Node "between" (10,52) [], Node "repR" (52,61) []])
						:((Node "iRep" (2,45) [Node "repL" (1,10) [], Node "between" (10,35) [], Node "repR" (35,44) []]):[])))
						
invertedRepeatsM :: SignalM Feature
invertedRepeatsM = liftF $ invertedRepeats

tandemRepeats :: String -> String -> SignalM Feature
tandemRepeats m name = do { (Node s (start,end) parts) <- (mergeM $ tandem1 $ (motifM m "")); return (Node name (start,end) parts)}

longerThan :: Int -> (Feature -> Bool)
longerThan l = (\(Node s (start,end) _) -> (end-start > l))

isFlanking :: Feature -> Bool
isFlanking (Node _ (_,_) parts) = if (length parts) < 2 then False else 
									let (Node _ (_,_) parts1) = parts!!1
									in if (length parts1) < 3 then False else
										let 
											(Node _ (start1,end1) _) = parts1!!0
											(Node _ (start2,end2) _) = parts1!!2
										in (end1-start1 < 5) && (end2-start2 < 5)
										
start = motifM "ATG" "ATG"
stop = (motifM "TAG" "TAG") `mplus` (motifM "TAA" "TAA") `mplus` (motifM "TGA" "TGA")
tpaseDomain = allMotifsM "CGTTATCTA" "tpase"
tpase = start `followedBy` tpaseDomain `followedBy` stop
tpases = mergeMany1 tpase

inner1 = tpases `followedBy` (tandemRepeats "GC" "GCreps")
inner2 = (invertedRepeatsM,2) `havingB` inner1
inner3 = (directRepeatsM,2) `havingB` inner2

dnaTransposon = sat inner3 isFlanking

main :: IO ()
main = do 
	-- args <- getArgs
	-- let x = dna (args !! 0)
	let x = (dna "GGAGGCATGCGGCGGCGGCGATAATCGGCGAATGCGGATGCGGCGTTATCTAGGTTAGGACTGCGCGCCGTTTACGGCGCGAGCTGGGCTTTGATTCGCCGCCGCCGTCCTCGCCGATTATCGCCGCCGCCTTTCTGCT")
	let y = (dna "GGATGTTATGTTTATGGG")
	let z = (dna "GGATGCGCGCGCTATGTTTATGGG")
	putStrLn $ prettyShow $ getValues $ search dnaTransposon x (1,-1)
	
