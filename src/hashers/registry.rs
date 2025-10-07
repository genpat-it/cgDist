// registry.rs - Hasher registry for managing available hashers

use std::collections::HashMap;
use super::traits::AlleleHasher;
use super::{Crc32Hasher, Sha256Hasher, Md5Hasher, SequenceHasher, HammingHasher};

/// Registry for available hashers - Production Ready
pub struct HasherRegistry {
    hashers: HashMap<String, Box<dyn AlleleHasher>>,
}

impl HasherRegistry {
    pub fn new() -> Self {
        let mut registry = Self {
            hashers: HashMap::new(),
        };
        
        // Register built-in hashers
        registry.register_hasher("crc32", Box::new(Crc32Hasher));
        registry.register_hasher("sha256", Box::new(Sha256Hasher));
        registry.register_hasher("md5", Box::new(Md5Hasher));
        registry.register_hasher("sequence", Box::new(SequenceHasher));
        registry.register_hasher("hamming", Box::new(HammingHasher));
        
        registry
    }
    
    /// Register a new hasher
    pub fn register_hasher(&mut self, name: &str, hasher: Box<dyn AlleleHasher>) {
        self.hashers.insert(name.to_string(), hasher);
    }
    
    /// Get a hasher by name
    pub fn get_hasher(&self, name: &str) -> Option<&dyn AlleleHasher> {
        self.hashers.get(name).map(|h| h.as_ref())
    }
    
    /// Check if a hasher exists
    pub fn has_hasher(&self, name: &str) -> bool {
        self.hashers.contains_key(name)
    }
    
    /// List all available hashers
    pub fn list_hashers(&self) -> Vec<(&str, &str)> {
        self.hashers.values()
            .map(|h| (h.name(), h.description()))
            .collect()
    }
    
    /// Get all hasher names
    pub fn get_hasher_names(&self) -> Vec<&str> {
        self.hashers.keys().map(|s| s.as_str()).collect()
    }
}

impl Default for HasherRegistry {
    fn default() -> Self {
        Self::new()
    }
}