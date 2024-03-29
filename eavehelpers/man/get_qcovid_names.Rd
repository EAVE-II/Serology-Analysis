% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/serology_analysis.R
\name{get_qcovid_names}
\alias{get_qcovid_names}
\title{Helper function that will extract the QCOVID names from a dataframe
after applying filtering on a condition,
and requiring there to be >n events with condition}
\usage{
get_qcovid_names(df, var = insufficient_response, condition = "Yes", n = 5)
}
\arguments{
\item{df}{input dataframe object}

\item{var}{the name of the variable in the input dataframe to be checked}

\item{condition}{the value of the var to be filtered on (default 'Yes')}

\item{n}{the number of events required to be present for the filtered condition  (default 5)}
}
\value{
a list of QCOVID names
}
\description{
It is used to find a QCOVID variables that have >n events after applying a filtering
}
