function  [signal, signalSE]  = testDataSize(num)


signal = zeros(num,120001);

for k=1:num
    signalSE(k).signal = complex(rand(15,120001),rand(15,120001));
end


if labindex == 1

else
    signal = [];
    signalSE = [];
end
