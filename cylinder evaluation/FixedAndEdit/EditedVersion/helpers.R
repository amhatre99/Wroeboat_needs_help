centered.seq = function(from, to, by)
{
  n = ceiling((to-from) / by)
  return (sapply(1:(n+2), function(i) {(to+from)/2-by*(n/2+0.5)+by*(i-1)}))
}

full.seq = function(from, to, by)
{
  b = seq(from, to, by)
  if (b[length(b)] < to) {
    b = c(b, b[length(b)]+by)
  }
  return (b)
}