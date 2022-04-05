% generate test matrices, similar 2D poisson

ns = 3;
nb = 2;

nt = ns;
n = nt*ns + nb;
diag = 2*ones(ns,1);
off_diag = -1*ones(ns,1);

A_ns = spdiags([off_diag, diag, off_diag], [-1,0,1], ns, ns);
%full(A)
I = speye(ns);

A = kron(A_ns, I) + kron(I, A_ns);
%full(A)
B = -rand(nb, ns*nt);
D = -rand(nb, nb) + (n+1)*speye(nb,nb);

Q = [A, B'; B, D];
%full(Q)

%eig(full(Q))

% generate rhs
n_rhs = 1;

rhs = rand(n, nb);

X = Q\rhs;

res_norm = norm(Q*X-rhs)

