include("./forward_test/main_body.jl");

reference=JSWAP.readmat("./forward_test/v3_1_reference.mat","data");
pass_test_forward=norm(v3-reference,2)/norm(reference,2)<.0001;
