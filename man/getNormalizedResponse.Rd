\name{getNormalizedResponse}
\alias{getNormalizedResponse}
\title{
	It combines Retention time Correction Scores (RCS) and Grouping Scores (GS)
}
\description{
	This function does unity based normalization on Retention time Correction 
	Scores (RCS) as well as Grouping Scores (GS). 
}
\usage{
getNormalizedResponse(response)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
	A List of all responses calculated by getRGTVValues for all experiments of an
	Design of Experiment
}
}
\details{
	Grouping Score (GS)	is calculated by:\cr
	'good groups'^2/'bad groups\cr\cr
	
	For all RCS and GS values unitiy based normalization is done. For every 
	experiment within the DoE these two values are added together and returned.
}
\value{
	A vector with RTGV values 
}
\author{
Gunnar Libiseller
}
\note{
	Since RCS and GS can be within completely different ranges, normalization has
	to be done to prevent an excessive influence of either RCS or GS.
}

\seealso{
  \code{\link{getRGTVValues}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (response) 
{
    good_groups <- sapply(response, "[[", "good_groups")
    bad_groups <- sapply(response, "[[", "bad_groups")
    bad_groups_bool <- bad_groups == 0
    bad_groups[bad_groups_bool] <- 1
    good_groups[bad_groups_bool] <- good_groups[bad_groups_bool] + 
        1
    group_ratio <- good_groups^2/bad_groups
    ARTS <- 1/sapply(response, "[[", "mean_rel_rt_diff")
    ARTS_penalty <- 1/sapply(response, "[[", "retcor_done")
    ARTS <- ARTS/ARTS_penalty
    norm_group_ratio <- (group_ratio - min(group_ratio))/(max(group_ratio) - 
        min(group_ratio))
    norm_ARTS <- (ARTS - min(ARTS))/(max(ARTS) - min(ARTS))
    return(norm_group_ratio + norm_ARTS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{method}
