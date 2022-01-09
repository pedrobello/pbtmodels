#-----------------Tested Functions - Mac OS and Windows
#
#' A function that calculates the required time to a desired fraction in cumulative curves and respective rate.
#'
#' This function allows you to calculate the time to a desired cumulative fraction and respective germination rate. Use this function on raw data to avoid loss of points closer to the desired cumulative fraction.
#' @param data time course and cumulative dataset. Several treatments can be used at once as long as it respects the template and column names provided. A column with time in hours (CumTime) + a column with cumulative fractions (CumFraction) are required with at least one additional column for relevant treatment (e.g., germination temperature or water potential)
#' @param fraction from 0 to 1 used to calculate the time required for that level to be obtained in the cumulative time course. Standard value is 0.5 (50 percent), to calculate the time to 50 percent germination (T50) and respective germination rate (GR50). fraction level can be entered and be used for calculation and change column name.
#' @param treat_1,treat_2,treat_3,treat_4,treat_5 are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here. These column names do not need to be informed in case the provided template file is used to organize the data.
#' @keywords Tx, GRx, germination speed, germination rate
#' @importFrom dplyr group_by_at
#' @importFrom dplyr mutate
#' @importFrom dplyr tally
#' @importFrom magrittr %>%
#' @export
#' @examples
#' calcspeed(Mydata)
calcspeed <- function(data, fraction, treat_1, treat_2, treat_3, treat_4, treat_5) {
  #Define the informed treatments (columns) to be grouped when calling the function or set standard to all columns besides CumFraction and CumTime.
  if (missing(treat_5)) { #treat_5 not informed
    if (missing(treat_4)) { #treat_4 not informed
      if (missing(treat_3)) { #treat_3 not informed
        if (missing(treat_2)) { #treat_2 not informed
          if (missing(treat_1)) { #treat_1 not informed
            treatcolnames <- c("treat.ID", "treat.desc", "treat.aging.time", "treat.priming.wp", "treat.priming.temp", "treat.priming.duration", "Germ.wp", "Germ.temp", "Germ.promoter.dosage", "Germ.inhibitor.dosage")
          } else {treatcolnames <- c(treat_1)}
        } else {treatcolnames <- c(treat_1, treat_2)}
      } else {treatcolnames <- c(treat_1, treat_2, treat_3)}
    } else {treatcolnames <- c(treat_1, treat_2, treat_3, treat_4)}
  } else {treatcolnames <- c(treat_1, treat_2, treat_3, treat_4, treat_5)}

  #Define the informed fraction for time and rate calculation. Standard is 0.5 -> 50 percent. Also populates new column name base in the fraction used.
  if (missing(fraction)) { #fraction not informed
    frac <- 0.5
    fracspeed_lbl <- "T50"
    fracrate_lbl <- "GR50"
  } else {
    frac <- fraction
    fracspeed_lbl <- paste("T", (frac * 100), sep = "")
    fracrate_lbl <- paste("GR", (frac * 100), sep = "")
  }

  # Calculate the time to the fraction selected using linear interpolation (approx function) after it calculates the inverse of time to generate the rate.
  treatments <- data %>% dplyr::group_by_at(treatcolnames) %>%
    dplyr::mutate(Tx = approx(CumFraction, CumTime, xout = frac, ties = "ordered")$y,
                  GRx = 1 / approx(CumFraction, CumTime, xout = frac, ties = "ordered")$y)

  treatcolnames <- c(treatcolnames, "Tx", "GRx")

  # Separate all treatments without germination time courses
  treatments <- treatments %>% dplyr::group_by_at(treatcolnames) %>% dplyr::tally()

  #Replaces the column name with the used fraction instead.
  names(treatments)[names(treatments) == "Tx"] <- fracspeed_lbl
  names(treatments)[names(treatments) == "GRx"] <- fracrate_lbl

  return(treatments)
}



#' A Function to clean cumulative curves on dataset.
#'
#' This function removes repetitive cumulative fractions/percentages, keeping only the initial presence of the value
#' @param data object with the raw cumulative data that needs to be removed.
#' @param treat_1,treat_2,treat_3,treat_4,treat_5 are the names of the treatment columns to separate the dataset. The time course cumulative curves will be grouped for each distinct treatment that should be informed here.
#' @keywords Clean cumulative fraction repetitive percentage
#' @importFrom dplyr distinct
#' @export
#' @examples cleandata(mydata,"treat.desc")
#' cleandata(mydata,"treat.desc")
cleandata <- function(data, treat_1, treat_2, treat_3, treat_4, treat_5) {
  #Clean Repetitive Percentages (keeps only initial presence of value)

  if (missing(data)) { #data object not informed
    print("Informe the data object.")
  } else {
    treatdata <- data
    if (missing(treat_1)) { #treatment 1 not informed
      print("Informed treatment for factor.")
    } else {
      t_1 <-  treat_1
      if (missing(treat_2)) { #treatment 2 not informed
        t_2 <- ""
      } else {
        t_2 <- paste(",", treat_2, sep = "")
      }
      if (missing(treat_3)) { #treatment 3 not informed
        t_3 <- ""
      } else {
        t_3 <- paste(",", treat_3, sep = "")
      }
      if (missing(treat_4)) { #treatment 2 not informed
        t_4 <- ""
      } else {
        t_4 <- paste(",", treat_4, sep = "")
      }
      if (missing(treat_5)) { #treatment 2 not informed
        t_5 <- ""
      } else {
        t_5 <- paste(",", treat_5, sep = "")
      }

      treatdataclean <- eval(parse(text = paste("dplyr::distinct(treatdata,", t_1, t_2, t_3, t_4, t_5, ", CumFraction, .keep_all = TRUE)", sep = "")))

      return(treatdataclean)
    }
  }
}

#----------------------New Development - Under Testing
