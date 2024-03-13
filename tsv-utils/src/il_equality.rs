pub trait ILEquality {
    fn switch_lj(self) -> Self;
}

impl ILEquality for String {
    fn switch_lj(self) -> Self {
        self.chars()
            .map(|character| match character {
                'L' => 'J',
                'J' => 'L',
                _ => character,
            })
            .collect()
    }
}

#[inline]
pub fn il_equality(char1: u8, char2: u8, i_and_l_equal: bool) -> bool {
    println!("{char1}, {char2}");
    if i_and_l_equal { // L is replaced by J in our custom alphabet
        char1 == char2 || ((char1 == b'I' && char2 == b'J') || (char1 == b'J' && char2 == b'I'))
    } else {
        char1 == char2
    }
}
