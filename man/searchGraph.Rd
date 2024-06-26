% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph.R
\name{searchGraph}
\alias{searchGraph}
\title{Perform Depth-first search from all nodes in a graph}
\usage{
searchGraph(graph, depthLimit = 20)
}
\arguments{
\item{graph}{the graph to search.}
}
\value{
a list of "traversals", which is a list of `Edge` objects in order they were traversed, for all nodes. See \link[chromolyse]{depthFirstSearch}.
}
\description{
Perform Depth-first search from all nodes in a graph
}
