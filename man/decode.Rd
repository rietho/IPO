\name{decode}
\alias{encode}
\alias{decode}
\alias{decodeAll}
\title{
	Decodes values from ranges of -1 to 1 to ranges specified.
}
\description{
	Encode and decode values that are in a range of -1 to 1 into a  
	specified range.
}
\usage{
	encode(value, bounds)
	decode(value, bounds)
    decodeAll(values, params)
}

\arguments{
  \item{value}{
	A value 
}
  \item{values}{
	A vector with values in the range [-1,1]
}
  \item{bounds}{
	A vector of two values defining the lower and upper bound of a range.
}
  \item{params}{
	A list where evere list-item consist of two values defining a lower
	and an upper bound.
}
}
\details{
	Decodes a values from ranges of -1 to 1 to ranges specified.

      	A function used to decode values that are in a range of -1 to 1 into a  
	specified range. For every value a list item with lower and upper bound
	has to be supplied.

	A function used to encode values that are in a specified range into a  
	range between -1 to 1.

      }
\value{
  decode: The encoded value.
  decodeAll:	A vector of decoded values.
}
\author{
	Gunnar Libiseller
}
\examples{

	decode(0, c(10, 20))
	decode(-0.5, c(10, 20))
	decode(1, c(10, 20))
	
	bounds <- c(10, 20)
	encode(decode(1, bounds), bounds)

        ## Multiple values:
	values <- c(-1, -0.25, 0, 0.75)
	params <- getDefaultXcmsSetStartingParams()
    type_params <- typeCastParams(params)
	
	decodeAll(values, type_params$to_optimize)

        ## Combination of encode and decode
	encode(15, c(10, 20))
	encode(10, c(10, 20))
	encode(5, c(1, 5))

	bounds <- c(1,5)
	decode(encode(5, bounds), bounds)

}
\keyword{method}
