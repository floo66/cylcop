.onLoad <- function(libname, pkgname) {
  if(.Platform$OS.type=="windows")
    sound::setWavPlayer(system.file("extdata", "wv_player.exe", package = "cylcop")) else
      sound::setWavPlayer("afplay")
}
