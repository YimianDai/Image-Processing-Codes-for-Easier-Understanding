function psnr = imgpsnr(img,nimg)

img  = double(img);
nimg = double(nimg);
mse  = sum(sum((nimg-img).^2)) / prod(size(img));
psnr = 20 * log10(255 / sqrt(mse));
% psnr = 20 * log10(max(max(img)) / sqrt(mse));