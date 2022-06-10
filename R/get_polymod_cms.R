library(socialmixr)



get_polymod_cms = function(breaks = c(c(0),seq(9,69, 10), c(120))){
  
  
  breaks = replace(breaks, breaks==Inf, 120)
  
  
  data("polymod")
  
  
  cms_pmd = contact_matrix(polymod, age.limits = breaks , countries = 'GB', n = 100)
  
  
  mr <- Reduce("+", lapply(cms_pmd$matrices, function(x) {x$matrix})) / length(cms_pmd$matrices)
  
  msd <- Reduce("+", lapply(cms_pmd$matrices, function(x) {abs(mr - x$matrix)})) / length(cms_pmd$matrices)
  
  cmd_pmd = setNames(melt(mr), c('Var1', 'Var2', 'cms'))
  cmd_pmd_sd = setNames(melt(msd), c('Var1', 'Var2', 'csd'))
  cmd_pmd = data.table::data.table(merge(cmd_pmd, cmd_pmd_sd, by=c('Var1', 'Var2')))
  
  cmd_pmd_all = data.table()
  
  for(srnd in 1:100){
    cmd_pmd_all = rbind(cmd_pmd_all, cmd_pmd[, sr:=srnd])
  }
  
  cmd_pmd_all
}


