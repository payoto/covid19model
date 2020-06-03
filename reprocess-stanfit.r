library(optparse)
library(data.table)
source("utils/log-and-process.r")


# Commandline options and parsing
parser <- OptionParser()
parser <- add_option(parser, c("-D", "--directory"), default=NULL,
                        help="Load all fits in directory X")
parser <- add_option(parser, c("-a","--all"), action="store_true",
                        help="Loads all fits in directory './'")
parser <- add_option(parser, c("-p","--pattern"), default="",
                        help="Finds files which match the specified pattern.")

cmdoptions <- parse_args(parser, args = commandArgs(trailingOnly = TRUE), 
    positional_arguments = TRUE)

if(!is.null(cmdoptions$options$all)) {
    cmdoptions$options$directory = "./"
}
if(is.null(cmdoptions$options$pattern)) {
    cmdoptions$options$pattern = ""
}

if(length(cmdoptions$args) > 0) {
    stanfit_files = cmdoptions$args
    if (!is.null(cmdoptions$options$directory)){
        message("WARNING: Ignoring directory as files specified on command line.")
    }
} else if (!is.null(cmdoptions$options$directory) 
    || !is.null(cmdoptions$options$pattern)) {
    if (is.null(cmdoptions$options$directory)){
        cmdoptions$options$directory = './'
    }
    stanfit_files = list.files(
        path=cmdoptions$options$directory,
        pattern=paste0(".*", cmdoptions$options$pattern, ".*-stanfit.Rdata$"),
        recursive = TRUE,
        full.names = TRUE
    )
} else {
    message("WARNING: no file specified")
    stanfit_files <- c()
}

print(stanfit_files)
error_log = data.frame(
    "error"=rep("", length(stanfit_files)),
    "file"=stanfit_files,
    stringsAsFactors=F
)
i = 0
for (stanfit_file in stanfit_files) {
    i = i+ 1
    err_mark = tryCatch({
      process_stanfit_file(stanfit_file)
      
    }, error = function(e) {
      message(paste("Processing Error:  ",e))
      return(paste(e))
    })
    if (!is.null(err_mark)){
      error_log$error[i] = err_mark
    }
}
message("______________________________________________")
message("________________ERROR_LOG_____________________")
message("______________________________________________")
print(error_log)