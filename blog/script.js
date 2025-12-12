console.log("System initialized. Welcome, user.");

// Optional: Add simple typewriter effect if elements are waiting
document.addEventListener('DOMContentLoaded', () => {
    const subtitles = document.querySelectorAll('.subtitle');
    subtitles.forEach((el, index) => {
        el.style.opacity = 0;
        setTimeout(() => {
            el.style.opacity = 1;
        }, (index + 1) * 500);
    });
});

