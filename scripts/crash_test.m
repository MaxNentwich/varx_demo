
for i = 1:1000

    x = randn(4,1e3,1e3);
    y = randn(4,1e3);
    
    z = tensorprod(x,y,[1 3],[1 2])';

end