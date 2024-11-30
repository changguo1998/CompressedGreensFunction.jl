using CompressedGreensFunction

grid = [
    0  1  0  0  0;
    0  0  1  0  0;
    0  0  0  1 -1.;
]

edge2grid = [
    1 2;
    1 3;
    2 3;
    1 4;
    2 4;
    3 4;
    1 5;
    2 5;
    3 5;
] |> permutedims

face2edge = [
    1 2 3;
    1 4 5;
    2 4 6;
    3 5 6;
    1 7 8;
    2 7 9;
    3 8 9;
] |> permutedims

elem2face = [
    1 2 3 4;
    1 5 6 7;
] |> permutedims

m = CompressedGreensFunction.Mesh(grid, edge2grid, face2edge, elem2face)
