The goal of this project is to implement a domain specific language using two different metaprogramming strategies. Here I have implemented a domain specific language for certain biological applications. 

Biologists and computational biologists are frequently required to look for different types of signals in DNA, RNA and protein sequences. These signals vary from simple substring motif to complex structures like repeats, inverted repeats or certain probabilistic models like hidden markov models defining the signals. Various algorithms are used to search for these signals but combining these signals is often a huge pain in the programmer's part. So, it is useful to create a domain specific language that abstracts the creation of these signals in a modular way and provides useful combinators that allow a programmer to combine multiple signals to create new signals. 

Here, I start by designing a language that could express a signal like this:

```
Signal tpase = ** TPase Pfam protein domain ** 
Signal tir = ** Terminal Inverted Repeat of length 15 to 270 **
Signal tsd = ** Target site duplication of TWA **
Signal pif = ** PIF Pfam protein domain **
Signal gcRich = ** GC-rich region **
Signal harbinger = (tsd left)
                    :> ((tir left) `hasSignal` gcRich)
                    :> tpase
                    :> (zeroOrExactlyOne pif)
                    :> (tir right)
                    :> (tsd right)
```
(:>) stands for "followed by". 

This is a real example of a specification for a type of biological element called a DNA transposon (sub-class "Harbinger"). The signal specifications can involve basic (or complex) string algorithms like KMP search, palindrome search, dynamic programming (local, global alignment), etc., or some machine learning algorithms like Hidden Markov Models, etc. In the example above, the Pfam protein domains are specified using "Profile Hidden Markov Models" (PHMM). 

