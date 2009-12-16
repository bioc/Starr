\name{match_ac}
\alias{match_ac}
\title{Exact String matching using the Aho-Corasick algorithm}
\description{
The Algorithm is designed to match a given set of sequences to a text. This function is used by remap for mapping reporter sequences to the genome.
}
\usage{
match_ac(dictionary, text, complementary=FALSE, reverse=FALSE, reverse_complementary=FALSE)
}
\arguments{
\item{ dictionary}{Set of sequences, given as a character vector}
\item{text}{At least one text, to match the sequences to.}
\item{complementary}{If TRUE, the sequences are searched in the complementary strand of the text}
\item{reverse}{If TRUE, the sequences are searched in the reverse strand of the text}
\item{reverse_complementary}{If TRUE, the sequences are searched in the reverse complementary strand of the text}
}
\author{ Benedikt Zacher \email{zacher@lmb.uni-muenchen.de}}
\keyword{internal}
