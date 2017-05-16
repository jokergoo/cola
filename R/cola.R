

#' A bottle of cola
#'
#' @export
#' @importFrom  crayon red green yellow blue magenta cyan silver
#' @importFrom  crayon bgRed bgGreen bgYellow bgBlue bgMagenta bgCyan bgWhite
#' @importFrom  crayon reset
cola = function() {
	cola_logo = 
"        _        
      .!.!.      
       ! !       
       ; :       
      ;   :      
     ;_____:     
     ! cola!     
     !_____!     
     :     :     
     :     ;     
     .'   '.     
     :     :     
      '''''      
"
    bg_color = sample(c("bgRed", "bgGreen", "bgYellow", "bgBlue", "bgMagenta", "bgCyan", "bgWhite"), 1)
	text_color = sample(c("red", "green", "yellow", "blue", "magenta", "cyan", "silver"), 1)

    fun = getFromNamespace(text_color, ns = "crayon")
    bg_fun = getFromNamespace(bg_color, ns = "crayon")
    
    lines = readLines(textConnection(cola_logo))
    lines = lines[-length(lines)]
    lines2 = NULL
    add_bg = runif(1) > 0.5
    for(i in seq_along(lines)) {
        if(add_bg) {
            lines2 = c(lines2, bg_fun(fun(lines[i])))
        } else {
            lines2 = c(lines2, fun(lines[i]))
        }
        lines2 = c(lines2, reset("\n"))    
    }

    cat(paste(lines2, collapse = ""))
}
