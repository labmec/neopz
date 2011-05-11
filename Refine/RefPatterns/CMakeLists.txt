
file(GLOB refpatterns *.rpt)
file(GLOB text *.txt)
install(FILES ${refpatterns} ${text} DESTINATION include/refpatterns)

