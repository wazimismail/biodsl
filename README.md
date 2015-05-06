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
`(:>)` stands for "followed by". 

This is a real example of a specification for a type of biological element called a DNA transposon (sub-class "Harbinger"). The signal specifications can involve basic (or complex) string algorithms like KMP search, palindrome search, dynamic programming (local, global alignment), etc., or some machine learning algorithms like Hidden Markov Models, etc. In the example above, the Pfam protein domains are specified using "Profile Hidden Markov Models" (PHMM). 

I have implemented the DSL using a naive approach (**biodsl.hs**) and tried implementing it using a modular monadic approach (**mbiodsl.hs**). 

- I have defined a data structure called `Feature`, which is nothing but a tree (of intervals). This data structure will be used to represent the structure and sub-structure of the result obtained by applying a `DNA` to a `Signal`. 

- `Signal` is a small DSL with a couple of combinators:
  1. the "followed-by" combinator (`:>`), which searches for the left signal followed by the right signal and builds the output `Feature` containing three parts (left, between and right). 
  2. the "contains" combinator (`:#:`), which searches for an outer signal whose part (specified by an integer) contains an inner signal, and builds the sub-structure of the output `Feature` in multiple parts. 

- I have implemented an interpreter 
  `findFeature :: Signal -> DNA -> Feature`. 
  This interpreter has the implementation for the two combinators, which involves de-structuring the output of one `Signal` and then constructing the combined result of the new `Signal`. 

- I have also implemented a bunch of signals to work with, which include `motif`, `invertedRepeats` and `gcRich`. 
- The rest of the code is helper functions and main-function for testing. 

In the second strategy, I have implemented this naive DSL using monadic modular approach (**mbiodsl.hs**). By defining the output type of `findFeature` as an `Identity` monad, wrapped by an `ErrorT` transformer, it becomes easy to handle the cases when the signals are not found and when this information should be propagated through the combinators. It reduces a lot of pattern matching under the hood. 

I would like to take this project forward by implementing the combinators within the monadic framework, to improve the functionality and modularity. But due to time constraints, I have presented only till this implementation. 

Advantages of modular monadic approach:
- It is easy to plug-in a certain type of language feature or side-effect using monad transformers, when they are available. 
- The flow of information is intuitive and clear. 

Disadvantages of modular monadic approach:
- Implementing some new monad transformers could involve a lot of plumbing under the hood. 

