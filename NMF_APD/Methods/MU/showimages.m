X = importdata('componentsAll30.arff', ','); 
for k=1:30
    img=zeros([28 28]);
    for i=1:28
        img(i,:) = X(k,(i-1)*28 + 1:i*28);  
    end
    subplot(6,5,k), subimage(img);
end;