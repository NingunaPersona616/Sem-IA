clear all
close all
clc

img = imread("display_ref.jpg");

[~, ~, P] = readBarcode(img, "QR-CODE");
P

img = insertShape(img, "FilledCircle", [P(1, :) 20], "Color", "red", "Opacity", 1);
img = insertShape(img, "FilledCircle", [P(2, :) 20], "Color", "blue", "Opacity", 1);
img = insertShape(img, "FilledCircle", [P(3, :) 20], "Color", "green", "Opacity", 1);

imshow(img)