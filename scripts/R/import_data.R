library(data.table)

f_get_parameter <- function(name_file, name_para) {
  # remove extension
  name_file <- sub("\\.csv$", "", name_file)
  
  # split only when '-' is followed by a letter (keeps negatives like -2.0)
  list_parameter <- unlist(strsplit(name_file, "-(?=[a-zA-Z])", perl = TRUE))[-1]
  
  # extract matching parameter
  parameter <- list_parameter[startsWith(list_parameter, paste0(name_para, "="))]
  if (length(parameter) == 0) {
    return(list(name_var = name_para, value_var = NA_character_))
  }
  
  name_var  <- strsplit(parameter, "=")[[1]][1]
  value_var <- strsplit(parameter, "=")[[1]][2]
  
  if (grepl("^\\[.*\\]$", value_var)) {
    # remove square brackets, split on comma, trim spaces
    value_list <- trimws(unlist(strsplit(gsub("^\\[|\\]$", "", value_var), ",")))
    name_list  <- paste0(name_var, seq_along(value_list))
    return(list(name_var = name_list, value_var = value_list))
  } else {
    return(list(name_var = name_var, value_var = trimws(value_var)))
  }
}

parse_suffix_number <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(tolower(as.character(x)))
  m <- regexec("^([+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+))(k?)$", x)
  mt <- regmatches(x, m)[[1]]
  if (length(mt) == 0) stop("Input not in recognised format.")
  number <- as.numeric(mt[2])
  mult <- if (mt[3] == "k") 1e3 else 1
  number * mult
}

try_parse_value <- function(x) {
  x <- trimws(as.character(x))
  xl <- tolower(x)
  if (xl %in% c("true","false")) return(as.logical(xl))
  if (xl %in% c("na","nan")) return(NA)
  num <- suppressWarnings(as.numeric(x))
  if (!is.na(num)) return(num)
  maybe <- tryCatch(parse_suffix_number(x), error = function(e) NA_real_)
  if (!is.na(maybe)) return(maybe)
  x
}

f_import_data <- function(name_file, listVar = c(), parse_value = TRUE, as_data_frame = FALSE) {
  ls <- list.files(pattern = "\\.csv$")

  if (length(name_file) != 0) {
    for (i in name_file) {
      ls <- grep(pattern = i, ls, value = TRUE, fixed = TRUE)
    }
  }
  
  if (!length(ls)) stop("No files matched your filters.")
  
  data <- data.table()
  
  for (file_name in ls) {
    print(file_name)
    
    dataTemp <- fread(file = file_name)
    
    for (var in listVar) {
      param_info <- f_get_parameter(file_name, var)
      name_vars  <- param_info$name_var
      value_vars <- param_info$value_var
      
      if (all(is.na(name_vars)) || all(is.na(value_vars))) next
      
      for (i in seq_along(name_vars)) {
        raw <- value_vars[i]
        val <- if (parse_value) try_parse_value(raw) else as.character(raw)
        dataTemp[, (name_vars[i]) := rep(val, .N)]
      }
    }
    
    data <- rbindlist(list(data, dataTemp), fill = TRUE, use.names = TRUE)
  }
  if (as_data_frame) return(as.data.frame(data))
  data
}