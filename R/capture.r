capture <- function(z,n){
	aft <- bef <- rep(1,2^n-1)
	aa <- bb <- m <- rep(1,n)
	for(i in (1-n):-1){
		m[n+i] <- z$fit[1]+z$fit[2^(n+i-1)+1]
		aft <- aft*(2-rep(rep(1:2,rep(2^(n+i-1),2)),2^(-i))[1:(2^n-1)])
		aa[n+i] <- sum(aft*z$fit[1:(2^n-1)])
		bef <- bef*(2-rep(rep(1:2,rep(2^(-i),2)),2^(n+i-1))[1:(2^n-1)])
		bb[1-i] <- sum(bef*z$fit[1:(2^n-1)])}
	aa[n] <- bb[1] <- z$fit[1]
	cc <- aa[c(n,1:(n-1))]
	m[n] <- cc[n]
	ph <- z$fit[1]/m
	nest <- aa*m*bb/z$fit[1]^2
	phih <- aa*m/(cc*z$fit[1])
	if(max(ph[2:(n-1)])-min(ph[2:(n-1)])<=10^(-4)){
		nest <- nest*ph/ph[2]
		ph[1] <- ph[2]
		ph[n] <- ph[2]
		phih[n] <- z$fit[1]/(m[n]*ph[n])}
	phih[1] <- 1
	phih[n] <- phih[n]+(m[n-1]/m[n]-1)*10^(-4)*(phih[n-1]<=1-10^(-4))
	bh <- nest-phih*nest[c(n,1:(n-1))]
	bh[1] <- nest[1]
	zz <- cbind(1:n,nest,phih,ph,bh)
	colnames(zz) <- c("i","N(i)","Phi(i-1)","P(i)","B(i-1)")
	zz}

setup <- expression({
	p1 <- as.numeric(gl(2,1,2^n))-1
	pbd <- p2 <- as.numeric(gl(2,2,2^n))-1
	p3 <- as.numeric(gl(2,4,2^n))-1
	b2 <- i1 <- p1*p2
	d1 <- i2 <- p2*p3
	pb <- pbd+p3
	if(n>=4){
		p4 <- as.numeric(gl(2,8,2^n))-1
		d2 <- i3 <- p3*p4
		b3 <- b2*p3
		d1 <- d2*p2
		pbd <- pbd+p3
		pb <- pbd+p4}
	if(n>=5){
		p5 <- as.numeric(gl(2,16,2^n))-1
		d3 <- i4 <- p4*p5
		b4 <- b3*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p4
		pb <- pbd+p5}
	if(n>=6){
		p6 <- as.numeric(gl(2,32,2^n))-1
		d4 <- i5 <- p5*p6
		b5 <- b4*p5
		d3 <- d4*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p5
		pb <- pbd+p6}
	if(n>=7){
		p7 <- as.numeric(gl(2,64,2^n))-1
		d5 <- i6 <- p6*p7
		b6 <- b5*p6
		d4 <- d5*p5
		d3 <- d4*p4
		d2 <- d3*p3
		d1 <- d2*p2
		pbd <- pbd+p6
		pb <- pbd+p7}
	pw <- rep(1,2^n)
	pw[2^n] <- 0
	pd <- pbd+p1
	pc <- pb+p1})

