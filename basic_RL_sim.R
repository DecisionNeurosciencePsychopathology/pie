basic_RL_sim=function(in_values,rew,seg_chosen,alpha) {
  
  out_values=in_values
  out_values=in_values[seg_chosen]+alpha*(rew-in_values[seg_chosen])
  
  return(out_values)
}