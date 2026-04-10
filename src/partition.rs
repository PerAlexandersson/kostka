//! Parsing helpers for partition and weight inputs.

use combinatoric_core::Partition;

/// Parse a comma-separated list of non-negative integers into a Partition.
pub fn parse_partition(s: &str) -> Result<Partition, String> {
    let mut parts: Vec<u32> = s.split(',')
        .map(|x| x.trim().parse::<u32>())
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("Invalid partition: {}", e))?;
    parts.sort_unstable_by(|a, b| b.cmp(a));
    Ok(Partition::new(parts))
}

/// Parse a comma-separated list into a weight vector (composition); zeros allowed, order kept.
pub fn parse_weight(s: &str) -> Result<Vec<u32>, String> {
    s.split(',')
        .map(|x| {
            x.trim()
                .parse::<u32>()
                .map_err(|e| format!("Invalid weight: {}", e))
        })
        .collect()
}
