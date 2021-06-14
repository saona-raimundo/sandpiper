fn main() {
    let args: Vec<String> = std::env::args().collect();
    let start: usize = args[1].parse().unwrap();

    println!("Hola!");
    println!("Starting from: {}", start);
}
