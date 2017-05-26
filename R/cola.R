
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

# == title
# A bottle of cola
# 
# == details
# Simply provide you a bottle of cola.
#
# The ASCII art is from http://ascii.co.uk/art/coke
#
# == example
# cola()
cola = function() {

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
        lines2 = c(lines2, crayon::reset("\n"))    
    }

    cat(paste(lines2, collapse = ""))
}
