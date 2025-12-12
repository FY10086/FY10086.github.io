document.addEventListener('DOMContentLoaded', () => {
    console.log("MANBO SYSTEM INITIALIZED");

    const startBtn = document.getElementById('start-btn');
    const startScreen = document.getElementById('start-screen');
    const audio = document.getElementById('bgm');

    startBtn.addEventListener('click', () => {
        // Play Audio
        audio.play().then(() => {
            console.log("Audio playing");
        }).catch(err => {
            console.error("Audio playback failed:", err);
        });

        // Hide Start Screen with a transition
        startScreen.style.opacity = '0';
        startScreen.style.transition = 'opacity 1s ease';
        setTimeout(() => {
            startScreen.style.display = 'none';
        }, 1000);

        // Start Chaos
        startChaos();
    });

    function startChaos() {
        const symbols = ['曼波', '🐟', '🐈', '哈基米', 'MAMBO', 'Huh?', '🤔', '🔥', '💊'];
        const colors = ['#ff00ff', '#00ffff', '#ffff00', '#ff0000', '#00ff00'];

        function createFloater() {
            const el = document.createElement('div');
            el.classList.add('floater');
            el.innerText = symbols[Math.floor(Math.random() * symbols.length)];
            
            // Random position
            el.style.left = Math.random() * 90 + 5 + 'vw';
            
            // Random size
            el.style.fontSize = (Math.random() * 3 + 1) + 'rem';
            
            // Random color
            el.style.color = colors[Math.floor(Math.random() * colors.length)];
            
            // Random duration
            el.style.animationDuration = (Math.random() * 3 + 2) + 's';

            document.body.appendChild(el);

            // Remove after animation
            setTimeout(() => {
                el.remove();
            }, 5000);
        }

        // Create a floater every 300ms
        setInterval(createFloater, 300);
    }

    // Glitch effect on header click (always active)
    const header = document.querySelector('h1');
    header.addEventListener('click', () => {
        document.body.style.filter = 'invert(1) hue-rotate(180deg)';
        setTimeout(() => {
            document.body.style.filter = 'none';
        }, 100);
    });
});
