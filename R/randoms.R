#' posteriorSummary function
#' I still need to document all of these things properly
#' @export
posteriorSummary <- function( paramSampleVec){
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  return( list( mean =  meanParam , median =  medianParam , mode = modeParam ))
}

#' @export
FigCount<-function(){

  thisEnv <- environment()

  me<-list(

    thisEnv = thisEnv,

    getEnv = function(){
      return(get("thisEnv",thisEnv))
    },

    init = function(){
      assign("figlist",list(),thisEnv)
      return(assign("figs",0,thisEnv))
    },

    newfig = function(label,caption){
      figs = get("figs",thisEnv)
      figs = figs + 1
      assign("figs",figs,thisEnv)
      figlist = get("figlist",thisEnv)
      figlist[[label]]$caption = caption
      figlist[[label]]$number = figs
      return(assign("figlist",figlist,thisEnv))
    },

    ref = function(label){
      figlist = get("figlist",thisEnv)
      return(figlist[[label]]$number)
    },

    cap = function(label){
      require(glue)
      figlist = get("figlist",thisEnv)
      return(glue("Figure {figlist[[label]]$number}: {figlist[[label]]$caption}"))
    }

  )

  assign("this",me,envir = thisEnv)

  class(me) <- append(class(me),"FigCount")
  return(me)
}

#' @export
TabCount<-function(){

  thisEnv <- environment()

  me<-list(

    thisEnv = thisEnv,

    getEnv = function(){
      return(get("thisEnv",thisEnv))
    },

    init = function(){
      assign("tablist",list(),thisEnv)
      return(assign("tabs",0,thisEnv))
    },

    newtab = function(label,caption){
      tabs = get("tabs",thisEnv)
      tabs = tabs + 1
      assign("tabs",tabs,thisEnv)
      tablist = get("tablist",thisEnv)
      tablist[[label]]$caption = caption
      tablist[[label]]$number = tabs
      return(assign("tablist",tablist,thisEnv))
    },

    ref = function(label){
      tablist = get("tablist",thisEnv)
      return(tablist[[label]]$number)
    },

    cap = function(label){
      require(glue)
      tablist = get("tablist",thisEnv)
      return(glue("Table {tablist[[label]]$number}: {tablist[[label]]$caption}"))
    }

  )

  assign("this",me,envir = thisEnv)

  class(me) <- append(class(me),"TabCount")
  return(me)
}

#' @export
'%!in%' <- function(x,y)!('%in%'(x,y))

#' @export
zScore.by <- function(data,index_var,measure_var){
  index = dplyr::enexpr(index_var)
  measure = dplyr::enexpr(measure_var)

  data %>% group_by(!!index) %>% summarise(sd = sd(!!measure), mean = mean(!!measure)) -> z.tmp
  merge.data.frame(data,z.tmp) %>% mutate(z = (!!measure - mean)/sd)
}

#' @export
send.alert<-function(message){
  require(glue)
  path = paste(system.file(package = "rara"),"emailalert.py", sep = "/")
  vals = read.csv("~/.mail",header = F, sep = "")
  email = as.character(vals[1,])
  password = as.character(vals[2,])
  command <- glue('python {path} "{message}" "{password}" "{email}"')
  system(command)
}
