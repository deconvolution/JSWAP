@parallel function Dx_1(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dx_1(in);
    return nothing;
end

@parallel function Dy_1(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dy_1(in);
    return nothing;
end

@parallel function Dz_1(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dz_1(in);
    return nothing;
end
