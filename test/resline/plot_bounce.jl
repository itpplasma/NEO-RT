using PyPlot

data = readdlm("bounce_zeroorder.out")

rc("text", usetex=true)
rc("figure", figsize=[2.8, 2.2])
rc("lines", linewidth=0.5)
rc("lines", markersize=3)

plot(sqrt(data[:,2]).*cos(data[:,4]), sqrt(data[:,2]).*sin(data[:,4]),"k:")
xlabel("\$x\$")
ylabel("\$y\$")
title("\$y\$")
tight_layout()
#savefig("test.png")

# a = 3
# b = "bla"
#
# ccall((:__common_MOD_disp, "libdriftorbit.so"), Int32, (Ptr{Cstring},Ptr{Cfloat}),
#        &b,&a)
