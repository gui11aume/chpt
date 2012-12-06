chpt <- function(x, n_CPU="auto", verbose=TRUE) {

   # Validate input type or stop with meaningful error message.
   if (!is.list(x)) {
      if(!is.vector(x)) {
         stop("'x' must be a vector or a list of vectors")
      }
      else {
         x <- list(x)
      }
   }
   if (!n_CPU=="auto" && !is.numeric(n_CPU)) {
      stop("'n_CPU' must be \"auto\" or a number")
   }
   if (!all(sapply(x, is.vector))) {
      stop("all the elements of 'x' must be vectors")
   }
   ref_length = length(x[[1]])
   for (this_length in lapply(x, length)) {
      if (any(this_length != ref_length)) {
         stop("all the matrices in 'x' must have same length")
      }
   }
   # Make sure values in 'x' are double.
   for (i in 1:length(x)) {
      storage.mode(x[[i]]) <- "double"
   }

   # Assign automatic variables and coerce to proper type.
   n_CPU <- as.integer(ifelse(n_CPU == "auto", 0, n_CPU))
   verbose <- as.logical(verbose)

   return(.Call("chpt_R_call", x, n_CPU, verbose))

}
