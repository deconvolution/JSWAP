include("./forward_test/main_body.jl");

#reference=JSWAP.readmat("./forward_test/v3_1_reference.mat","data");
#pass_test_forward=norm(v3-reference,2)/norm(reference,2)<.0001;
tt=(v3[20,20,20]-(-2.3666761533221374e-15))/(-2.3666761533221374e-15)<10.0^-5;
tt2=(v3[30,30,20]-(-3.442374644804789e-6))/(-3.442374644804789e-6)<10.0^-5;
tt3=(v3[30,33,30]-(0.15315952419043086))/(0.15315952419043086)<10.0^-5;
pass_test_forward=tt && tt2 && tt3;
