function V = position_mod_LE(V,box_size,strain)

Vx = V(:,1);
Vy = V(:,2);

ind = Vy >= box_size;
Vx(ind) = mod(Vx(ind) - strain*box_size,box_size);

ind = Vy >= 0 & Vy < box_size;
Vx(ind) = mod(Vx(ind),box_size);

ind = Vy < 0;
Vx(ind) = mod(Vx(ind) + strain*box_size,box_size);


Vy = mod(Vy,box_size);

V = [Vx,Vy];


end