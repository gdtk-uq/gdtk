import std.math;
import std.range;
import std.random;
import std.algorithm;
import plt = matplotlibd.pyplot;

void main() {
    simple();
    color();
    polar();
    subplots();
}

void bad() {
    plt.plot("hoge");
    plt.show();
}

void simple() {
    auto x = iota(0, 2.05, 0.05).map!(x => x * PI);
    auto y = x.map!(sin);
    
    plt.plot(x, y, "r-", ["label": "$y=sin(x)$"]);
    plt.xlim(0, 2 * PI);
    plt.ylim(-1, 1);
    plt.legend();
    plt.savefig("simple.png");
    plt.clear();
}

void color() {
    const n = 100;
    auto x = iota(n);
    auto y = x[];
    double[n][n] z;

    foreach (i; 0..n)
        foreach (j; 0..n)
            z[i][j] = i + j;
    plt.contourf(x, y, z, 64, ["cmap": "hsv"]);
    plt.colorbar();
    plt.savefig("color.png");
    plt.clear();
}

void subplots() {
    auto x = iota(0, 2 * PI + 0.05, 0.05);

    plt.subplot(221);
    plt.plot(x, x.map!(sin));
    plt.xlim(0, 2 * PI);
    plt.ylim(-1, 1);

    plt.subplot(222);
    plt.plot(x, x.map!(cos));
    plt.xlim(0, 2 * PI);
    plt.ylim(-1, 1);

    plt.subplot(223);
    plt.plot(x, x.map!(i => sin(i) * exp(-0.4 * i)));
    plt.xlim(0, 2 * PI);
    plt.ylim(-1, 1);

    plt.subplot(224);
    plt.plot(x, x.map!(i => cos(i) * exp(-0.4 * i)));
    plt.xlim(0, 2 * PI);
    plt.ylim(-1, 1);

    plt.savefig("subplots.png");
    plt.clear();
}

void polar() {
    
    auto r = iota(0, 1.001, 0.001);
    auto theta = r.map!(i => i * 32 * PI);
    auto area = r.map!(i => i * 2000);
    
    plt.subplot(111, ["projection": "polar"]);
    plt.scatter(theta, r, ["c": r], ["s": area], ["cmap": "hsv"], ["alpha": 0.25]);
    plt.savefig("polar.png");
    plt.clear();
}

