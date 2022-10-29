# if(.Platform$OS.type=="windows")
#   sound::setWavPlayer(system.file("extdata", "wv_player.exe", package = "cylcop")) else
#     sound::setWavPlayer("afplay")
cylcop.env <- new.env(parent = emptyenv())
assign("silent", TRUE, envir=cylcop.env)
