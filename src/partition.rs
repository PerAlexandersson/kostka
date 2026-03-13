/// A partition stored as a Vec<u32> of parts in weakly decreasing order, no trailing zeros.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Partition(pub Vec<u32>);

impl Partition {
    pub fn new(mut parts: Vec<u32>) -> Self {
        while parts.last() == Some(&0) {
            parts.pop();
        }
        Partition(parts)
    }

    pub fn empty() -> Self {
        Partition(vec![])
    }

    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    pub fn num_parts(&self) -> usize {
        self.0.len()
    }

    pub fn size(&self) -> u32 {
        self.0.iter().sum()
    }

    /// Part at row r (0-indexed), 0 if r >= num_parts.
    pub fn part(&self, r: usize) -> u32 {
        self.0.get(r).copied().unwrap_or(0)
    }

    /// Returns true if self ⊆ other (containment of Young diagrams).
    pub fn contained_in(&self, other: &Partition) -> bool {
        self.0.iter().enumerate().all(|(i, &p)| p <= other.part(i))
    }

    /// Conjugate partition.
    pub fn conjugate(&self) -> Partition {
        if self.0.is_empty() {
            return Partition::empty();
        }
        let max_part = self.0[0] as usize;
        let mut conj = vec![0u32; max_part];
        for &p in &self.0 {
            for j in 0..p as usize {
                conj[j] += 1;
            }
        }
        Partition(conj)
    }

    /// All partitions of n with at most max_parts parts, each part <= max_part.
    pub fn all_of_size_bounded(n: u32, max_parts: usize, max_part: u32) -> Vec<Partition> {
        let mut result = Vec::new();
        Self::_generate(n, max_parts, max_part, &mut vec![], &mut result);
        result
    }

    fn _generate(
        remaining: u32,
        parts_left: usize,
        max_part: u32,
        current: &mut Vec<u32>,
        result: &mut Vec<Partition>,
    ) {
        if remaining == 0 {
            result.push(Partition(current.clone()));
            return;
        }
        if parts_left == 0 {
            return;
        }
        let last = current.last().copied().unwrap_or(max_part);
        let upper = last.min(remaining).min(max_part);
        for p in (1..=upper).rev() {
            current.push(p);
            Self::_generate(remaining - p, parts_left - 1, max_part, current, result);
            current.pop();
        }
    }

    /// All partitions of n.
    pub fn all_of_size(n: u32) -> Vec<Partition> {
        Self::all_of_size_bounded(n, n as usize, n)
    }

    pub fn display(&self) -> String {
        if self.0.is_empty() {
            return "∅".to_string();
        }
        self.0.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
    }
}

impl std::fmt::Display for Partition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display())
    }
}

/// Parse a comma-separated list of non-negative integers into a Partition.
pub fn parse_partition(s: &str) -> Result<Partition, String> {
    let parts: Result<Vec<u32>, _> = s.split(',')
        .map(|x| x.trim().parse::<u32>())
        .collect();
    parts
        .map(Partition::new)
        .map_err(|e| format!("Invalid partition: {}", e))
}

/// Parse a comma-separated list into a weight vector (composition); zeros allowed, order kept.
pub fn parse_weight(s: &str) -> Result<Vec<u32>, String> {
    s.split(',')
        .map(|x| x.trim().parse::<u32>().map_err(|e| format!("Invalid weight: {}", e)))
        .collect()
}
