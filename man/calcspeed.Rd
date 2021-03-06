% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicFunctions.R
\name{calcspeed}
\alias{calcspeed}
\title{A function that calculates the required time to a desired fraction in cumulative curves and respective rate.}
\usage{
calcspeed(data, fraction, treat_1, treat_2, treat_3, treat_4, treat_5)
}
\arguments{
\item{data}{time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)}

\item{fraction}{from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Standard value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). fraction level can be entered and be used for calculation and change column name.}

\item{treat_1, treat_2, treat_3, treat_4, treat_5}{are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.}
}
\description{
This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
}
\examples{
calcspeed(Mydata)
}
\keyword{GRx,}
\keyword{Tx,}
\keyword{germination}
\keyword{rate}
\keyword{speed,}
