source_gitfile <- function(filename){
  url = paste0("https://raw.githubusercontent.com/zhimeir/adaptive_knockoff_paper/master/",filename,".R")
  script = GET(url = url)
  script = content(script,"text")
  eval(parse(text = script),envir= .GlobalEnv)
}

file_vec <- c("R/generic_getorder","R/plot_vanilla","R/plot_ordering")
getfile <- sapply(file_vec,source_gitfile)


# run simulation 1 
ParamsRowIndex=1
#source simulation 1
res =filter_EM_getorder(W,z,alpha = 0.2,offset=1,mute = FALSE,df=2)
plot_ordering(res,nonzero = nonzero,start_index = 650)
ggsave(paste0("~/Dropbox/Adaptive Knockoff Paper/plots/adaordering1.pdf"),plot = last_plot(),width = 10,height = 4)
plot_vanilla(W,nonzero = nonzero,start_index = 650,alpha=0.2)
ggsave(paste0("~/Dropbox/Adaptive Knockoff Paper/plots/vordering1.pdf"),plot = last_plot(),width = 10,height = 4)

# run simulation 2
ParamsRowIndex=1
#source simulation 2
res =filter_EM_getorder(W,z,alpha = 0.2,offset=1,mute = FALSE,df=2)
plot_ordering(res,nonzero = nonzero,start_index = 1300)
ggsave(paste0("~/Dropbox/Adaptive Knockoff Paper/plots/adaordering2.pdf"),plot = last_plot(),width = 10,height = 4)
plot_vanilla(W,nonzero = nonzero,start_index = 1300,alpha=0.2)
ggsave(paste0("~/Dropbox/Adaptive Knockoff Paper/plots/vordering2.pdf"),plot = last_plot(),width = 10,height = 4)
