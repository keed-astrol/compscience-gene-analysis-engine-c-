import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import subprocess
import threading
import os
import sys
import queue

try:
    from PIL import Image, ImageTk
    print("PIL/Pillow is installed successfully")
except ImportError:
    print("PIL/Pillow is not installed. Please install it using: pip install pillow")
    sys.exit(1)

# Determine the executable name based on the platform
if sys.platform.startswith('win'):
    exe_name = "projectfinal.exe"
else:
    exe_name = "./projectfinal"

class BackgroundApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DNA Sequence Analysis")
        self.geometry("620x500")
        # Make window resizable
        self.resizable(True, True)
        self.minsize(620, 500)  # Set minimum size
        
        # Variables to hold file paths
        self.input_file = tk.StringVar()
        self.output_file = tk.StringVar()
        
        # Status queue for thread communication
        self.status_queue = queue.Queue()
        
        # Set a fallback background color
        self.configure(bg='navy')
        
        # Create a container frame that will expand with the window
        self.main_container = tk.Frame(self)
        self.main_container.pack(fill=tk.BOTH, expand=True)
        
        # Load and set background image
        self.load_background()
        
        # Create UI elements
        self.create_widgets()
        
        # Start checking the queue
        self.after(100, self.check_queue)
        
        # Bind resize event
        self.bind("<Configure>", self.on_resize)
    
    def load_background(self):
        image_path = r"C:\Users\ADMIN\source\repos\compscienceproject\IMAGE.jpg"
        
        try:
            if os.path.exists(image_path):
                print(f"Found image at: {image_path}")
                
                # Load the original image
                self.original_image = Image.open(image_path)
                print(f"Original image size: {self.original_image.size}")
                
                # Create the initial background
                self.update_background_image()
            else:
                print(f"Image not found at: {image_path}")
                self.prompt_for_image()
        except Exception as e:
            print(f"Error loading image: {str(e)}")
            self.use_fallback_background()
    
    def update_background_image(self):
        try:
            # Get current window size
            width = self.winfo_width()
            height = self.winfo_height()
            
            # Use at least the minimum size
            width = max(width, 620)
            height = max(height, 500)
            
            print(f"Resizing image to: {width}x{height}")
            
            # Resize the image to match the current window size
            resized_image = self.original_image.resize((width, height), Image.LANCZOS)
            
            # Convert to PhotoImage
            self.bg_photo = ImageTk.PhotoImage(resized_image)
            
            # Check if bg_label exists, update it or create it
            if hasattr(self, 'bg_label'):
                self.bg_label.configure(image=self.bg_photo)
                self.bg_label.image = self.bg_photo  # Keep a reference
            else:
                # Create a label with the background image
                self.bg_label = tk.Label(self.main_container, image=self.bg_photo)
                self.bg_label.place(x=0, y=0, relwidth=1, relheight=1)
                self.bg_label.image = self.bg_photo  # Keep a reference
                
            print("Background image updated")
        except Exception as e:
            print(f"Error updating background: {str(e)}")
    
    def on_resize(self, event):
        # Only respond to main window resize events, not child widgets
        if event.widget == self:
            # Add a small delay to avoid excessive updates
            self.after_cancel(self.after_id) if hasattr(self, 'after_id') else None
            self.after_id = self.after(100, self.update_background_image)
    
    def prompt_for_image(self):
        print("Prompting user to select an image")
        messagebox.showinfo("Image Not Found", "Background image not found. Please select an image file.")
        new_image_path = filedialog.askopenfilename(
            title="Select Background Image",
            filetypes=[("Image files", "*.jpg *.jpeg *.png")]
        )
        if new_image_path:
            print(f"User selected image: {new_image_path}")
            try:
                self.original_image = Image.open(new_image_path)
                self.update_background_image()
            except Exception as e:
                print(f"Failed to load user-selected image: {str(e)}")
                self.use_fallback_background()
        else:
            print("No image selected by user")
            self.use_fallback_background()
    
    def use_fallback_background(self):
        print("Using fallback solid background")
        # Use a solid color background instead
        self.main_container.configure(bg="navy")
    
    def create_widgets(self):
        # Create a transparent frame for the form that will be positioned relative to window size
        self.form_container = tk.Frame(self.main_container)
        self.form_container.place(relx=0.5, rely=0.5, anchor=tk.CENTER, relwidth=0.9, relheight=0.8)
        
        # Create the semi-transparent form frame
        form_frame = tk.Frame(self.form_container, bg='black')
        form_frame.place(relx=0.5, rely=0.5, anchor=tk.CENTER, relwidth=1, relheight=0.5)
        
        # Style for the labels
        label_style = {'bg': 'black', 'fg': 'yellow', 'font': ('Segoe UI', 10, 'bold')}
        
        # Input file row
        input_label = tk.Label(form_frame, text="Input FASTA File:", **label_style)
        input_label.grid(row=0, column=0, sticky="e", padx=10, pady=15)
        
        entry_input = ttk.Entry(form_frame, textvariable=self.input_file, width=50)
        entry_input.grid(row=0, column=1, padx=10, pady=15)
        
        btn_input = ttk.Button(form_frame, text="Browse", command=self.select_input_file)
        btn_input.grid(row=0, column=2, padx=10, pady=15)
        
        # Output file row
        output_label = tk.Label(form_frame, text="Output CSV File:", **label_style)
        output_label.grid(row=1, column=0, sticky="e", padx=10, pady=15)
        
        entry_output = ttk.Entry(form_frame, textvariable=self.output_file, width=50)
        entry_output.grid(row=1, column=1, padx=10, pady=15)
        
        btn_output = ttk.Button(form_frame, text="Browse", command=self.select_output_file)
        btn_output.grid(row=1, column=2, padx=10, pady=15)
        
        # Run button
        btn_run = ttk.Button(form_frame, text="▶ Run Analysis", command=self.run_analysis)
        btn_run.grid(row=2, column=1, pady=20)
        
        # Status label
        self.status_label = tk.Label(form_frame, text="", bg='black', fg='yellow', font=('Segoe UI', 10, 'bold'))
        self.status_label.grid(row=3, column=0, columnspan=3, pady=5)
        
        # Configure the grid columns to expand with the window
        form_frame.columnconfigure(1, weight=1)
    
    def select_input_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
        if file_path:
            self.input_file.set(file_path)
    
    def select_output_file(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if file_path:
            self.output_file.set(file_path)
    
    def run_analysis(self):
        if not self.input_file.get() or not self.output_file.get():
            messagebox.showerror("Error", "Please select both input and output files.")
            return
        
        exe_path = os.path.join(os.path.dirname(__file__), exe_name)
        if not os.path.exists(exe_path):
            messagebox.showerror("Error", "C++ executable not found at " + exe_path)
            return
        
        self.status_label.config(text="Running analysis...")
        thread = threading.Thread(target=self.run_cpp_program)
        thread.daemon = True
        thread.start()
    
    def run_cpp_program(self):
        try:
            exe_path = os.path.join(os.path.dirname(__file__), exe_name)
            result = subprocess.run(
                [exe_path, self.input_file.get(), self.output_file.get()], 
                capture_output=True, 
                text=True
            )
            if result.returncode == 0:
                self.status_queue.put("Analysis complete. Results saved to " + self.output_file.get())
            else:
                self.status_queue.put("Error: " + result.stderr)
        except Exception as e:
            self.status_queue.put("Exception: " + str(e))
    
    def check_queue(self):
        try:
            while True:
                message = self.status_queue.get_nowait()
                self.status_label.config(text=message)
        except queue.Empty:
            pass
        self.after(100, self.check_queue)

# Create and run the application
if __name__ == "__main__":
    app = BackgroundApp()
    # Give the window a chance to initialize before setting up the background
    app.update_idletasks()
    app.after(100, app.update_background_image)
    app.mainloop()