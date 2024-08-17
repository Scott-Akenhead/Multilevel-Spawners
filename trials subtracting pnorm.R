# trials with subtracting pnorm()
x=0:15
residence = 2
a0 = dnorm(x, 5, 2);
a =  pnorm(x, 5, 2)
b0 = dnorm(x, 5+ residence, 2)
b =  pnorm(x, 5+ residence, 2)
d0 = (a-b)/residence
d= cumsum(d0)
plot(x,a, type='l', ylab="Spawners", xlab="Time",
     sub='black: Entry, blue: Exit, red: Presence/Residence',
     main='Entry = N(5,2), Residence = 2')
abline(v=c(5,5+residence), h=0.5, lty='dotted')
lines(x,a0)
lines(x,b0, col="blue");lines(x,b, col="blue")
lines(x,d0, col="red", lwd=3)
lines(x,d, col="red", lwd=3)

cat('a \n',round(a,2),'\n'); cat('a0 \n',round(a0,2),'\n');
cat('b \n',round(b,2),'\n'); cat('b0 \n',round(b0,2),'\n');
cat('d \n',round(d,2),'\n'); cat('d0 \n',round(d0,2),'\n');

Presence = function(days, mu, sigma, residence){
    a = (pnorm(days,mu,sigma) - pnorm(days,mu+residence,sigma))/ residence
    return(a)
}

t= c(3,5,7)
r1 = Presence (t, 5, 2, residence)
cat('days and presence \n',t,'\n', round(r1,3),'\n')
  # days and presence
  #     3     5     7
  # 0.068 0.171 0.171  checks out!
