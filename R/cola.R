

#' A bottle of cola
#'
#' @export
#' @importFrom  crayon red green yellow blue magenta cyan silver
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
	text_color = c("red", "green", "yellow", "blue", "magenta", "cyan", "silver")
	cat(getFromNamespace(sample(text_color, 1), ns = "crayon")(cola_logo))
}
