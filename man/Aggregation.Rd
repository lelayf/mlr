% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Aggregation.R
\name{Aggregation}
\alias{Aggregation}
\title{Aggregation object.}
\description{
An aggregation method reduces the performance values of the test
(and possibly the training sets) to a single value.
To see all possible implemented aggregations look at \code{\link{aggregations}}.

The aggregation can access all relevant information of the result after resampling
and combine them into a single value. Though usually something very simple
like taking the mean of the test set performances is done.

Object members:
\describe{
\item{id [\code{character(1)}]}{Name of the aggregation method.}
\item{name [\code{character(1)}]}{Long name of the aggregation method.}
\item{properties [\code{character}]}{Properties of the aggregation.}
\item{fun [\code{function(task, perf.test, perf.train, measure, group, pred)}]}{Aggregation function.}
}
}
\seealso{
\code{\link{makeAggregation}}
}

