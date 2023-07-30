/// Perform a max k-selection (e.g. select the 50 largest items) on an
/// array, in place. This algorithm works by building a bounded min heap
/// in the array. The `k` highest elements in the slice will be stored in
/// first `k` elements in the slice, but they will be stored in min-heap order,
/// not in sorted order. Otherwise, this is equivalent to performing a partial
/// sort
pub fn bounded_min_heapify<T: Ord>(slice: &mut [T], k: usize) {
    if slice.len() <= k {
        return;
    }

    // Build a heap in place for the first `k` elements
    for i in (0..k / 2).rev() {
        sift_down(&mut slice[..k], i);
    }

    debug_assert!(check_heap(&slice[..k]));

    // Scan the vector, "inserting" items into the heap if they are larger
    // than the smallest item in the heap. This performs the `k` select operation
    for i in k..slice.len() {
        if slice[i] > slice[0] {
            slice.swap(i, 0);
            sift_down(&mut slice[..k], 0);
            debug_assert!(check_heap(&slice[..k]));
        }
    }
}

fn check_heap<T: Ord>(slice: &[T]) -> bool {
    for i in 1..slice.len() {
        let parent = (i - 1) / 2;
        if slice[parent] > slice[i] {
            return false;
        }
    }
    true
}

fn sift_down<T: Ord>(slice: &mut [T], mut index: usize) {
    while let Some(left) = slice.get(index * 2 + 1) {
        let mut smallest = index;
        if left < &slice[smallest] {
            smallest = index * 2 + 1;
        }

        if let Some(right) = slice.get(index * 2 + 2) {
            if right < &slice[smallest] {
                smallest = index * 2 + 2;
            }
        }

        if smallest != index {
            slice.swap(smallest, index);
            index = smallest;
        } else {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Debug;

    use quickcheck_macros::quickcheck;

    use super::bounded_min_heapify;
    use super::check_heap;

    fn check<T: Ord + Clone + Debug>(mut data: Vec<T>, k: usize) {
        let k = k.min(data.len());
        let mut cloned = data.clone();
        // Stable sort the data
        cloned.sort_by(|a, b| b.cmp(&a));

        bounded_min_heapify(&mut data, k);

        // Take the heap part, and sort it
        let top_k = &mut data[..k];

        // Check that heap property is maintained, or that k == length of the data
        assert!(check_heap(top_k) || k == cloned.len());

        top_k.sort_by(|a, b| b.cmp(&a));
        assert_eq!(top_k, &mut cloned[..k]);
    }

    #[quickcheck]
    fn run_quickcheck(data: Vec<i32>, k: usize) {
        check(data, k);
    }

    #[test]
    fn smoke() {
        let asc = (0..500).collect::<Vec<_>>();
        let desc = (0..500).rev().collect::<Vec<_>>();
        check(asc, 50);
        check(desc, 50);
    }
}
