
fprintf('***************************************\n');
fprintf('2D Phase Unwrapping Demo\n');
fprintf('Please select the demo:\n');
fprintf('(1) No noise  , no ignored region\n');
fprintf(' 2. With noise, no ignored region\n');
fprintf(' 3. No noise  , with ignored region\n');
fprintf(' 4. With noise, with ignored region\n');
while (1)
    user_input = input('Your selection (1-4): ', 's');
    user_input = strip(user_input);

    % if the user does not supply anything, select the default
    if strcmp(user_input, '')
        fprintf('Demo 1 is selected\n');
        user_input = '1';
    end

    if length(user_input) == 1 && sum(user_input == '1234') == 1
        break;
    else
        fprintf('Invalid input\n');
    end
end


[X, Y] = meshgrid(linspace(-1, 1, 512) * 5);
img = -(X.*X + Y.*Y);
fprintf('Image size: %dx%d pixels\n', size(img,1), size(img,2));

% add noise
if any(user_input == '24')
    img = img + randn(size(X)) * 0.5;
end

% add an ignored region
if any(user_input == '34')
    img(end/4:3*end/4,end/4:3*end/4) = nan;
end

% wrap the image
wimg = wrapTo2Pi(img);

tic;
unwrap_img = unwrap_phase(wimg);
toc;

subplot(221);
pcolor(img);
shading flat;
set(gca, 'ydir', 'reverse');
title('Original phase');

subplot(222);
pcolor(wimg);
shading flat;
set(gca, 'ydir', 'reverse');
title('Wrapped phase');

subplot(223);
pcolor(unwrap_img);
shading flat;
set(gca, 'ydir', 'reverse');
title('Unwrapped phase');

subplot(224);
pcolor(wrapTo2Pi(unwrap_img));
shading flat;
set(gca, 'ydir', 'reverse');
title('Rewrap of unwrapped phase');
