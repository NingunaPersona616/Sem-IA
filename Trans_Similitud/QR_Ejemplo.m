clear all
close all
clc

img = imread("ref_3.png");

[~, ~, P] = readBarcode(img, "QR-CODE");

img = insertShape(img, "FilledCircle", [P(1, :) 10], "Color", "red", "Opacity", 1);
img = insertShape(img, "FilledCircle", [P(2, :) 10], "Color", "blue", "Opacity", 1);
img = insertShape(img, "FilledCircle", [P(3, :) 10], "Color", "green", "Opacity", 1);

imshow(img)