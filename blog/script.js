document.addEventListener('DOMContentLoaded', () => {
    console.log("MANBO SYSTEM INITIALIZED");

    const symbols = ['曼波', '🐟', '🐈', '哈基米', 'MAMBO', 'Huh?'];
    const colors = ['#ff00ff', '#00ffff', '#ffff00', '#ff0000', '#00ff00'];

    function createFloater() {
        const el = document.createElement('div');
        el.classList.add('floater');
        el.innerText = symbols[Math.floor(Math.random() * symbols.length)];
        
        // Random position
        el.style.left = Math.random() * 90 + 5 + 'vw';
        
        // Random size
        el.style.fontSize = (Math.random() * 2 + 1) + 'rem';
        
        // Random color
        el.style.color = colors[Math.floor(Math.random() * colors.length)];
        
        // Random duration
        el.style.animationDuration = (Math.random() * 5 + 3) + 's';

        document.body.appendChild(el);

        // Remove after animation
        setTimeout(() => {
            el.remove();
        }, 8000);
    }

    // Create a floater every 500ms
    setInterval(createFloater, 500);

    // Glitch effect on header click
    const header = document.querySelector('h1');
    header.addEventListener('click', () => {
        document.body.style.filter = 'invert(1)';
        setTimeout(() => {
            document.body.style.filter = 'none';
        }, 100);
    });
});
