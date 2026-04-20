/// Union-Find with path halving and union by rank.
/// Tracks cluster sizes and the running maximum cluster size.
pub struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
    size: Vec<u32>,
    pub max_size: u32,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
            size: vec![1; n],
            max_size: 1,
        }
    }

    /// Find root with path halving (iterative, no recursion).
    pub fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] as usize != x {
            let grandparent = self.parent[self.parent[x] as usize];
            self.parent[x] = grandparent;
            x = grandparent as usize;
        }
        x
    }

    /// Union x and y by rank. Returns true if they were in different components.
    pub fn union(&mut self, x: usize, y: usize) -> bool {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return false;
        }
        let (root, child) = if self.rank[rx] >= self.rank[ry] {
            (rx, ry)
        } else {
            (ry, rx)
        };
        self.parent[child] = root as u32;
        if self.rank[root] == self.rank[child] {
            self.rank[root] += 1;
        }
        self.size[root] += self.size[child];
        if self.size[root] > self.max_size {
            self.max_size = self.size[root];
        }
        true
    }

    /// Size of the component containing x.
    pub fn size_of(&mut self, x: usize) -> u32 {
        let root = self.find(x);
        self.size[root]
    }

    /// Reset to n singletons without reallocating.
    #[allow(dead_code)]
    pub fn reset(&mut self, n: usize) {
        for i in 0..n {
            self.parent[i] = i as u32;
            self.rank[i] = 0;
            self.size[i] = 1;
        }
        self.max_size = 1;
    }
}
